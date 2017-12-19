#include <cassert>
#include <cmath>
#include <cfloat>
#include <fstream>

#include "FindNClusters.hh"
#include "fftjet/SimplePredicates.hh"

namespace {
    struct LogScaleCalc : public fftjet::Functor1<double, fftjet::Peak>
    {
        inline double operator()(const fftjet::Peak& p) const
        {
            return log(p.scale());
        }
    };
}

FindNClusters::FindNClusters(const unsigned injets, const double factor,
                             const char* filename)
    : integFactor(factor), corrector(0), njets(injets)
{
    assert(integFactor > 0.0);
    std::ifstream in(filename, std::ios_base::in |
                               std::ios_base::binary);
    if (in.is_open())
    {
        fftjet::LinearInterpolator2d* response = fftjet::LinearInterpolator2d::read(in);
        if (response)
            corrector = new fftjet::JetMagnitudeMapper2d<fftjet::Peak>(
                *response, new LogScaleCalc(), true,
                response->xMin(), response->xMax(), 1000,
                response->yMax(), 1000);
    }
}

FindNClusters::~FindNClusters()
{
    delete corrector;
}

int FindNClusters::run(const fftjet::AbsClusteringTree<fftjet::Peak,long>& tree,
                       std::vector<fftjet::Peak>* peaks, PatternRecoInfo* info)
{
    assert(peaks);
    assert(info);

    info->jetCount = njets;
    bool contiguous = tree.clusterCountLevels(
        info->jetCount, fftjet::Always<true,fftjet::Peak>(),
        &info->minLevel, &info->maxLevel);
    if (info->minLevel == 0 && info->maxLevel == 0)
    {
        // First, try to find one more peak. Maybe, we are in a rare
        // situation when the number of peaks jumps by 2 from one
        // level to the other
        ++info->jetCount;
        contiguous = tree.clusterCountLevels(
            info->jetCount, fftjet::Always<true,fftjet::Peak>(),
            &info->minLevel, &info->maxLevel);
        if (info->minLevel == 0 && info->maxLevel == 0)
        {
            // It seems we ask for too many peaks.
            // Now, try to reduce the number of peaks.
            for (info->jetCount = njets-1; info->jetCount > 0; --info->jetCount)
            {
                contiguous = tree.clusterCountLevels(
                    info->jetCount, fftjet::Always<true,fftjet::Peak>(),
                    &info->minLevel, &info->maxLevel);
                if (info->minLevel || info->maxLevel)
                    break;
            }
        }
    }
    if (info->jetCount == 0)
        // Can't find any clusters in the clustering tree
        // (or there are always too many)
        return 1;

    info->foundLevel = (info->minLevel + info->maxLevel)/2;
    if (!contiguous)
    {
        // We are not guaranteed that the middle level
        // contains the number of clusters we really want
        if (tree.nClusters(info->foundLevel) != info->jetCount)
            for (unsigned delta=1; true; ++delta)
            {
                if (info->foundLevel-delta >= info->minLevel)
                {
                    if (tree.nClusters(info->foundLevel-delta) == info->jetCount)
                    {
                        info->foundLevel = info->foundLevel-delta;
                        break;
                    }
                }
                else if (info->foundLevel+delta <= info->maxLevel)
                {
                    if (tree.nClusters(info->foundLevel+delta) == info->jetCount)
                    {
                        info->foundLevel = info->foundLevel+delta;
                        break;
                    }
                }
                else
                    // We should never ever get here
                    assert(!"Bug detected in the clusterCountLevels function");
            }
    }

    // Get the peak data from this level
    double scale;
    long dummy;
    tree.getLevelData(info->foundLevel, &scale, peaks, &dummy);

    // From each peak magnitude, make a first guess about the jet momentum
    const double norm = integFactor*scale*scale;
    for (unsigned i=0; i<info->jetCount; ++i)
    {
        const double ptGuess = (*corrector)(
            norm*(*peaks)[i].magnitude(), (*peaks)[i]);
        if (ptGuess > 1.0/FLT_MAX)
            (*peaks)[i].setRecoScale(1.0/ptGuess);
        else
            (*peaks)[i].setRecoScale(FLT_MAX);
    }

    return 0;
}
