#include "../counting/feature_counter.hpp"
#include "../counting/record_types.hpp"
#include "../support/algorithms.hpp"
#include "../support/streams.hpp"
#include "../support/time.hpp"
#include <bitset>
#include <unordered_map>

namespace nustar::star2711b::solo::detail {

constexpr std::uint32_t uncolored_node = std::numeric_limits<std::uint32_t>::max();

void collapseUMIwith1MMlowHalf(std::uint32_t* umiArr, std::uint32_t umiArrayStride,
                               std::uint32_t umiMaskLow, std::uint32_t nU0, std::uint32_t& nU1,
                               std::uint32_t& nU2, std::uint32_t& nC,
                               std::vector<std::array<std::uint32_t, 2>>& vC);
void graphDepthFirstSearch(std::uint32_t n, std::vector<std::vector<std::uint32_t>>& nodeEdges,
                           std::vector<std::uint32_t>& nodeColor);
std::uint32_t graphNumberOfConnectedComponents(std::uint32_t N,
                                               std::vector<std::array<std::uint32_t, 2>> V,
                                               std::vector<std::uint32_t>& nodeColor);

std::uint32_t FeatureCounter::correctUmisGraph(const std::uint32_t nU0, UmiCode* umiArr,
                                               const bool readInfoRec, const bool nUMIyes,
                                               std::unordered_map<UmiCode, UmiCode>& umiCorr) {
    std::uint32_t nU1 = nU0;
    std::uint32_t nU2 = nU0;
    std::uint32_t graphN = 0;                            // number of nodes
    std::vector<std::array<std::uint32_t, 2>> graphConn; // node connections
    std::vector<std::uint32_t> graphComp; // for each node (color) - connected component number

    for (UnsignedWord iu = 0; iu < nU0 * umiArrayStride; iu += umiArrayStride)
        umiArr[iu + 2] = uncolored_node; // marks no color for graph

    qsort(umiArr, nU0, umiArrayStride * sizeof(std::uint32_t), funCompareNumbers<std::uint32_t>);
    collapseUMIwith1MMlowHalf(umiArr, umiArrayStride, options.umiMaskLow, nU0, nU1, nU2, graphN,
                              graphConn);

    // exchange low and high half of UMIs, re-sort, and look for 1MM again
    for (std::uint32_t iu = 0; iu < umiArrayStride * nU0; iu += umiArrayStride) {
        options.umiSwapHalves(umiArr[iu]);
    };
    qsort(umiArr, nU0, umiArrayStride * sizeof(std::uint32_t), funCompareNumbers<std::uint32_t>);
    collapseUMIwith1MMlowHalf(umiArr, umiArrayStride, options.umiMaskLow, nU0, nU1, nU2, graphN,
                              graphConn);

    std::uint32_t nConnComp = graphNumberOfConnectedComponents(graphN, graphConn, graphComp);
    nU1 += nConnComp;

    if (readInfoRec) {
        for (std::uint32_t ii = 0; ii < graphComp.size();
             ii++) { // for non-conflicting colors, need to fill the colors correctly
            if (graphComp[ii] == (std::uint32_t)-1)
                graphComp[ii] = ii;
        };

        const std::uint32_t bitTopMask = ~(1 << 31);
        std::vector<std::array<std::uint32_t, 2>> umiBest(graphN, {0, 0});
        std::unordered_map<UmiCode, std::uint32_t> umiCorrColor;
        for (std::uint32_t iu = 0; iu < umiArrayStride * nU0; iu += umiArrayStride) {
            // switch low/high to recover original UMIs
            options.umiSwapHalves(umiArr[iu]); // halves were swapped, need to return back to UMIs
            // find best UMI (highest count) for each connected component
            if (umiArr[iu + 2] == uncolored_node)
                continue; // UMI is not corrected
            std::uint32_t color1 = graphComp[umiArr[iu + 2]];
            std::uint32_t count1 = umiArr[iu + 1] & bitTopMask;
            if (umiBest[color1][0] < count1) {
                umiBest[color1][0] = count1;
                umiBest[color1][1] = umiArr[iu];
            };
            umiCorrColor[umiArr[iu]] = color1;
        };

        for (std::uint32_t iu = 0; iu < umiArrayStride * nU0; iu += umiArrayStride) {
            auto umi = umiArr[iu + 0];
            if (umiCorrColor.count(umi) > 0)
                umiCorr[umi] = umiBest[umiCorrColor[umi]][1];
        };
    };

    if (nUMIyes) { // this is not needed
        return nU1;
    } else {
        return 0;
    };
};
void collapseUMIwith1MMlowHalf(std::uint32_t* umiArr, std::uint32_t umiArrayStride,
                               std::uint32_t umiMaskLow, std::uint32_t nU0, std::uint32_t& nU1,
                               std::uint32_t& nU2, std::uint32_t& nC,
                               std::vector<std::array<std::uint32_t, 2>>& vC) {
    const std::uint32_t bitTop = 1 << 31;
    const std::uint32_t bitTopMask = ~bitTop;

    for (std::uint32_t iu = 0; iu < umiArrayStride * nU0; iu += umiArrayStride) { // each UMI
        std::uint32_t iuu = iu + umiArrayStride;
        for (; iuu < umiArrayStride * nU0; iuu += umiArrayStride) { // compare to all UMIs down

            std::uint32_t uuXor = umiArr[iu] ^ umiArr[iuu];

            if (uuXor > umiMaskLow)
                break; // upper half is different

            if (uuXor >> (__builtin_ctz(uuXor) / 2) * 2 >
                3)        // shift by even number of trailing zeros
                continue; //>1MM

            // 1MM UMI

            // graph coloring
            if (umiArr[iu + 2] == uncolored_node && umiArr[iuu + 2] == uncolored_node) { // no color
                // new color
                umiArr[iu + 2] = nC;
                umiArr[iuu + 2] = nC;
                ++nC;
                nU1 -= 2; // subtract the duplicated UMIs
            } else if (umiArr[iu + 2] == uncolored_node) {
                umiArr[iu + 2] = umiArr[iuu + 2];
                --nU1; // subtract the duplicated UMIs
            } else if (umiArr[iuu + 2] == uncolored_node) {
                umiArr[iuu + 2] = umiArr[iu + 2];
                --nU1;                                   // subtract the duplicated UMIs
            } else {                                     // both color
                if (umiArr[iuu + 2] != umiArr[iu + 2]) { // color conflict
                    vC.push_back({umiArr[iu + 2], umiArr[iuu + 2]});
                    // vC.push_back({umiArr[iuu+2],umiArr[iu+2]});
                };
            };

            //"directional" collapse from the UMI-tools by Smith, Heger and Sudbery (Genome Research
            //2017).
            if ((umiArr[iuu + 1] & bitTop) == 0 &&
                (umiArr[iu + 1] & bitTopMask) >
                    (2 * (umiArr[iuu + 1] & bitTopMask) - 1)) { // iuu is duplicate of iu
                umiArr[iuu + 1] |= bitTop;
                --nU2; // subtract the duplicated UMIs
            } else if ((umiArr[iu + 1] & bitTop) == 0 &&
                       (umiArr[iuu + 1] & bitTopMask) >
                           (2 * (umiArr[iu + 1] & bitTopMask) - 1)) { // iu is duplicate of iuu
                umiArr[iu + 1] |= bitTop;
                --nU2; // subtract the duplicated UMIs
            };
        };
    };
};

void graphDepthFirstSearch(std::uint32_t n, std::vector<std::vector<std::uint32_t>>& nodeEdges,
                           std::vector<std::uint32_t>& nodeColor) {
    for (const auto& nn : nodeEdges[n]) {
        if (nodeColor[nn] == (std::uint32_t)-1) { // node not visited
            nodeColor[nn] = nodeColor[n];
            graphDepthFirstSearch(nn, nodeEdges, nodeColor);
        };
    };
};

std::uint32_t graphNumberOfConnectedComponents(
    std::uint32_t N, std::vector<std::array<std::uint32_t, 2>> V,
    std::vector<std::uint32_t>& nodeColor) { // find number of connected components
    // N=number of nodes
    // V=edges, list of connected nodes, each pair of nodes listed once
    // simple recursive DFS

    nodeColor.resize(
        N,
        (std::uint32_t)-1); // new color (connected component) for each node (each original color)

    if (V.size() == 0)
        return N;

    std::vector<std::vector<std::uint32_t>> nodeEdges(N);
    for (std::uint32_t ii = 0; ii < V.size(); ii++) {
        nodeEdges[V[ii][0]].push_back(V[ii][1]);
        nodeEdges[V[ii][1]].push_back(V[ii][0]);
    };

    std::uint32_t nConnComp = 0;
    for (std::uint32_t ii = 0; ii < N; ii++) {
        if (nodeEdges[ii].size() == 0) { // this node is not connected, no need to check. Save time
                                         // because this happens often
            ++nConnComp;
        } else if (nodeColor[ii] == (std::uint32_t)-1) { // node not visited
            ++nConnComp;
            nodeColor[ii] = ii;
            graphDepthFirstSearch(ii, nodeEdges, nodeColor);
        };
    };
    return nConnComp;
};

} // namespace nustar::star2711b::solo::detail
