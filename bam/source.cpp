#include "../third_party/htslib/htslib/sam.h"
#include "io.hpp"
#include <algorithm>
#include <stdexcept>
namespace nustar::star2711b {
namespace {
struct FileCloser {
    void operator()(samFile* p) const {
        if (p)
            sam_close(p);
    }
};
struct HeaderCloser {
    void operator()(bam_hdr_t* p) const { bam_hdr_destroy(p); }
};
struct RecordCloser {
    void operator()(bam1_t* p) const { bam_destroy1(p); }
};
class BamSource final : public IReadSource {
    InputSpec spec_;
    std::size_t group_ = 0;
    std::uint64_t ordinal_ = 0, batch_ = 0;
    std::unique_ptr<samFile, FileCloser> file_;
    std::unique_ptr<bam_hdr_t, HeaderCloser> header_;
    std::unique_ptr<bam1_t, RecordCloser> record_{bam_init1()};
    kstring_t text_buffer_{0, 0, nullptr};
    bool cancelled_ = false, eof_ = false;
    void open() {
        file_.reset(sam_open(spec_.groups[group_].mate_paths[0].c_str(), "r"));
        if (!file_)
            throw std::runtime_error("cannot open SAM/BAM input");
        header_.reset(sam_hdr_read(file_.get()));
        if (!header_)
            throw std::runtime_error("invalid SAM/BAM header");
    }
    bool read(RawMate& mate, unsigned& flag) {
        const int result = sam_read1(file_.get(), header_.get(), record_.get());
        if (result == -1)
            return false;
        if (result < 0)
            throw std::runtime_error("truncated or invalid SAM/BAM record");
        const auto& b = *record_;
        flag = b.core.flag;
        if (b.core.l_qseq > 1 && bam_get_qual(&b)[0] == 255)
            throw std::runtime_error(
                "SAM/BAM input QUAL is missing (STAR requires sequence-length qualities)");
        if (b.core.l_qseq == 0)
            throw std::runtime_error("SAM/BAM input record has no sequence");
        mate.header = "@" + std::string(bam_get_qname(&b));
        mate.sequence.resize(b.core.l_qseq);
        mate.quality.emplace(b.core.l_qseq, '!');
        const auto* seq = bam_get_seq(&b);
        const auto* qual = bam_get_qual(&b);
        for (int i = 0; i < b.core.l_qseq; ++i) {
            mate.sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)];
            (*mate.quality)[i] = qual[0] == 255 ? '*' : char(qual[i] + 33);
        }
        if (flag & BAM_FREVERSE) {
            std::reverse(mate.sequence.begin(), mate.sequence.end());
            for (auto& c : mate.sequence) {
                switch (c) {
                case 'A':
                    c = 'T';
                    break;
                case 'T':
                    c = 'A';
                    break;
                case 'C':
                    c = 'G';
                    break;
                case 'G':
                    c = 'C';
                    break;
                default:
                    break;
                }
            }
            std::reverse(mate.quality->begin(), mate.quality->end());
        }
        // Solo consumes string tags. HTSlib handles framing and numeric arrays;
        // format only the auxiliary tail when the caller actually needs it.
        if (spec_.capture_aux) {
            auto& line = text_buffer_;
            line.l = 0;
            if (sam_format1(header_.get(), &b, &line) < 0) {
                throw std::runtime_error("invalid SAM/BAM auxiliaries");
            }
            std::string_view text(line.s, line.l);
            std::size_t pos = 0;
            for (unsigned column = 0; column < 11; ++column) {
                pos = text.find('\t', pos);
                if (pos == text.npos)
                    break;
                ++pos;
            }
            if (pos != text.npos)
                mate.header_extra = std::string(text.substr(pos));
        }
        return true;
    }

  public:
    explicit BamSource(const InputSpec& spec) : spec_(spec) { open(); }
    ~BamSource() override { free(text_buffer_.s); }
    Result<std::optional<ReadBatch>> next(std::uint64_t target) override {
        try {
            ReadBatch batch{BatchId{batch_++}, {}};
            std::uint64_t bytes = 0;
            while (!cancelled_ && !eof_ && bytes < target) {
                RawMate first;
                unsigned flag = 0;
                if (!read(first, flag)) {
                    if (++group_ == spec_.groups.size()) {
                        eof_ = true;
                        break;
                    }
                    open();
                    continue;
                }
                RawFragment raw{};
                raw.id = {std::uint32_t(group_), ordinal_++};
                raw.layout = spec_.layout;
                raw.format = FastxFormat::Fastq;
                raw.legacy_read_filter = (flag & 0x800) ? 'Y' : 'N';
                if (spec_.layout == Layout::PE) {
                    RawMate second;
                    unsigned other = 0;
                    if (!read(second, other) || first.header != second.header ||
                        !(((flag & 0x40) && (other & 0x80)) || ((flag & 0x80) && (other & 0x40))))
                        throw std::runtime_error(
                            "paired SAM/BAM requires adjacent records with the same QNAME and "
                            "complementary mate flags; collate input first");
                    // STAR loads mate 2 last into its per-fragment filter flag,
                    // regardless of the order of records in the input file.
                    const auto mate2_flag = (flag & 0x80) ? flag : other;
                    raw.legacy_read_filter = (mate2_flag & 0x800) ? 'Y' : 'N';
                    if (flag & 0x80)
                        std::swap(first, second);
                    raw.mates[1] = std::move(second);
                }
                raw.mates[0] = std::move(first);
                for (const auto& mate : raw.mates)
                    if (mate)
                        bytes += mate->header.size() + 2 * mate->sequence.size() +
                                 mate->header_extra.size();
                batch.fragments.push_back(std::move(raw));
            }
            if (batch.fragments.empty())
                return std::optional<ReadBatch>{};
            return std::optional<ReadBatch>{std::move(batch)};
        } catch (const std::exception& e) {
            return Failure{ErrorKind::Input, Stage::Input, {}, 1, e.what()};
        }
    }
    void cancel() noexcept override { cancelled_ = true; }
};
class Backend final : public IReadSourceBackend {
    Result<std::unique_ptr<IReadSource>> open(const InputSpec& spec) const override {
        try {
            return std::unique_ptr<IReadSource>(std::make_unique<BamSource>(spec));
        } catch (const std::exception& e) {
            return Failure{ErrorKind::Input, Stage::Input, {}, 1, e.what()};
        }
    }
};
} // namespace
const IReadSourceBackend& bamSourceBackend() {
    static Backend backend;
    return backend;
}
} // namespace nustar::star2711b
