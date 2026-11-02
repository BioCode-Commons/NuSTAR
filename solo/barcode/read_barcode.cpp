#include "read_barcode.hpp"
#include "../support/streams.hpp"

namespace nustar::star2711b::solo::detail {

ReadBarcode::ReadBarcode(Context& context_) : context_(context_), options(context_.options) {
    if (options.type == 0)
        return;

    for (std::uint32_t ii = 0; ii < stats.nStats; ii++)
        stats.V[ii] = 0;

    if (options.cbWLyes) {
        cbReadCountExact.resize(options.cbWLsize, 0);
    };

    for (std::uint32_t jj = 0; jj < 4; jj++) {
        homoPolymer[jj] = 0;
        for (std::uint32_t ii = 0; ii < options.umiL; ii++) {
            homoPolymer[jj] = (homoPolymer[jj] << 2) + jj;
        };
    };

    qualHist.fill(0);
};

void ReadBarcode::addCounts(const ReadBarcode& rfIn) {
    if (options.cbWLyes) {
        for (std::uint32_t ii = 0; ii < options.cbWLsize; ii++) {
            cbReadCountExact[ii] += rfIn.cbReadCountExact[ii];
        };
    };

    for (std::uint32_t ii = 0; ii < qualHist.size(); ii++)
        qualHist[ii] += rfIn.qualHist[ii];
};

void ReadBarcode::addStats(const ReadBarcode& rfIn) {
    for (std::uint32_t ii = 0; ii < stats.nStats; ii++)
        stats.V[ii] += rfIn.stats.V[ii];
};

void ReadBarcode::statsOut(std::ofstream& streamOut) {
    // streamOut << setw(50) << "CELL BARCODES IN READS:\n"
    for (std::uint32_t ii = 0; ii < stats.nStats; ii++) {
        streamOut << std::setw(50) << stats.names.at(ii) << std::setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};

} // namespace nustar::star2711b::solo::detail
