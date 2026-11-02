#pragma once

namespace nustar::star2711b::solo::detail {

inline constexpr char solo_defaults[] = R"defaults(
runMode                         alignReads
runThreadN                      1
readQualityScoreBase        33
soloType                    None
soloCBtype                  Sequence
soloCBwhitelist             -
soloCBstart                 1
soloCBlen                   16
soloUMIstart                17
soloUMIlen                  10
soloBarcodeReadLength       1
soloBarcodeMate             0
soloCBposition              -
soloUMIposition             -
soloAdapterSequence         -
soloAdapterMismatchesNmax   1
soloCBmatchWLtype           1MM_multi
soloInputSAMattrBarcodeSeq  -
soloInputSAMattrBarcodeQual  -
soloStrand                  Forward
soloFeatures                Gene
soloMultiMappers            Unique
soloUMIdedup                1MM_All
soloUMIfiltering            -
soloOutFileNames            Solo.out/          features.tsv barcodes.tsv        matrix.mtx
soloCellFilter              CellRanger2.2 3000 0.99 10
soloOutFormatFeaturesGeneField3    "Gene Expression"
soloCellReadStats           None
soloClusterCBfile           -
)defaults";

} // namespace nustar::star2711b::solo::detail
