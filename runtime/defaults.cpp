#include "defaults.hpp"

namespace nustar::star2711b {

StitchCommonConfig defaultStarStitchCommonConfig() {
    StitchCommonConfig c{};
    c.alignIntronMin = 21;
    c.alignIntronMax = 0;
    c.alignSJoverhangMin = 5;
    c.alignSJDBoverhangMin = 3;
    c.alignSJstitchMismatchNmax = {0, -1, 0, 0};
    c.alignSplicedMateMapLmin = 0;
    c.alignSplicedMateMapLminOverLmate = 0.66;
    c.alignTranscriptsPerWindowNmax = 100;
    c.alignTranscriptsPerReadNmax = 10000;
    c.alignEndsType = "Local";
    c.alignSoftClipAtReferenceEnds = true;
    c.insertion_flush_right = false;
    c.outFilterIntronMotifs = "None";
    c.outFilterIntronStrands = "RemoveInconsistentStrands";
    c.outSAMstrandField = "None";
    c.outFilterMismatchNoverLmax = 0.3;
    c.outFilterMultimapScoreRange = 1;
    return c;
}

SingleStitchConfig defaultStarSingleStitchConfig() {
    return SingleStitchConfig{defaultStarStitchCommonConfig(), {false, false}};
}

PairedStitchConfig defaultStarPairedStitchConfig() {
    PairedStitchConfig result{};
    result.common = defaultStarStitchCommonConfig();
    result.force_end_extension = {{{false, false}, {false, false}}};
    result.alignMatesGapMax = 0;
    result.protrusion_bases_max = 0;
    result.protrusion_concordant_pair = false;
    return result;
}

ScoreConfig defaultStarScoreConfig() {
    ScoreConfig result{};
    result.scoreGap = 0;
    result.scoreGapNoncan = -8;
    result.scoreGapGCAG = -4;
    result.scoreGapATAC = -8;
    result.scoreDelOpen = -2;
    result.scoreDelBase = -2;
    result.scoreInsOpen = -2;
    result.scoreInsBase = -2;
    result.scoreStitchSJshift = 1;
    result.sjdbScore = 2;
    result.scoreGenomicLengthLog2scale = -0.25;
    return result;
}


SeedConfig defaultStarSeedConfig() { return {50, 0, 1.0, 10000, 1000, 12, 5, 10, 0}; }
SelectionConfig defaultStarSelectionConfig() { return {1, 10, {}, "Old_2.4", "OneBestScore"}; }
ReadFilterConfig defaultStarReadFilterConfig() { return {0, 0.66, 0, 10, 0.66, 0.3, 1.0}; }
SamConfig defaultStarSamConfig() {
    SamConfig sam{};
    sam.enabled = true;
    sam.outSAMmode = "Full";
    sam.outSAMattributes = {"Standard"};
    sam.outSAMattrIHstart = 1;
    sam.outSAMmapqUnique = 255;
    sam.outSAMflagAND = 65535;
    sam.unmapped_within = true;
    return sam;
}

} // namespace nustar::star2711b
