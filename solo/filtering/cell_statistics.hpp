#pragma once
#include "../support/types.hpp"

namespace nustar::star2711b::solo::detail {

class CellStatistics {
public:
    std::array<UnsignedWord, 13> C;

    UnsignedWord& nCells = C[0];
    UnsignedWord& nReadInCells = C[1];
    UnsignedWord& medianReadPerCellUnique = C[2];
    UnsignedWord& meanReadPerCellUnique = C[3];
    UnsignedWord& nUMIinCells = C[4];
    UnsignedWord& medianUMIperCell = C[5];
    UnsignedWord& meanUMIperCell = C[6];
    UnsignedWord& nGeneInCells = C[7];
    UnsignedWord& medianGenePerCell = C[8];
    UnsignedWord& meanGenePerCell = C[9];
    UnsignedWord& nGeneDetected = C[10];
    UnsignedWord& nCellsSimple = C[11];
    UnsignedWord& nReadInCellsUnique = C[12];

    std::vector<std::uint32_t> nReadPerCellUnique, nGenePerCell;
    std::vector<bool> filtVecBool;

    void reset(UnsignedWord nCells) {
        C = {0};
        nReadPerCellUnique.clear();
        nReadPerCellUnique.reserve(16384);
        nGenePerCell.clear();
        nGenePerCell.reserve(16384);
        filtVecBool.clear();
        filtVecBool.resize(nCells, false);
    };
};

} // namespace nustar::star2711b::solo::detail
