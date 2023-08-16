#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <omp.h>

const int32_t IMPL = 1;

struct MaxSum
{
    int32_t maxSum;
    int32_t startIdx;
    int32_t endIdx;

    MaxSum(): maxSum(-1), startIdx(-1), endIdx(-1) {}
    MaxSum(int32_t _maxSum, int32_t _startIdx, int32_t _endIdx): maxSum(maxSum), startIdx(_startIdx), endIdx(_endIdx)
    {
    }

    //
    // Required for the reduction step in the parallelization.
    bool operator<(MaxSum& other) const
    {
        // Required when we are reducing results in Parallel Scan and
        // there are two threads that found the sequence with the maximal number of candies,
        // and one is unitialized, so we distinguish between one that finished and one that did not.
        if (maxSum == other.maxSum)
        {
            return startIdx < other.startIdx && endIdx < other.endIdx;
        }
        return maxSum < other.maxSum;
    }
};

uint32_t getNumberOfThreads(void)
{
    uint32_t numThreads;
    #pragma omp parallel
    {
    #pragma omp single
        numThreads = omp_get_num_threads();
    }
    return numThreads;
}

void readFile(int32_t& nCandies, std::vector<int32_t>& vHomeCandies)
{
    const std::string fileInputName = "input.txt";
    int32_t nHomes;

    std::ifstream fInput;
    fInput.open(fileInputName);
    if (!fInput.is_open())
    {
        std::cerr << "Wrong file name" << std::endl;
    }

    fInput >> nHomes;
    fInput >> nCandies;
    vHomeCandies.resize(nHomes);

    for (uint32_t idx = 0; idx < nHomes; idx ++)
    {
        fInput >> vHomeCandies[idx];
    }
    fInput.close();
}

// TODO: Compute prefix sums using parallelization.
void computePrefixSums(const std::vector<int32_t>& vHomeCandies, std::vector<int32_t>& vPrefSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    vPrefSum.resize(nHomes+1);
    int32_t curSum = 0;
    vPrefSum[0] = 0;
    for (uint32_t idx = 0; idx < nHomes; idx++)
    {
        curSum += vHomeCandies[idx];
        vPrefSum[idx+1] = curSum;
    }
}

void computeMaxSequenceSequential(const std::vector<int32_t>& vHomeCandies, const std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                   MaxSum& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    maxSum = MaxSum();

    for (int32_t idx = 1; idx <= nHomes; idx ++)
    {
        int32_t diff = std::max(vPrefSum[idx] - nCandies, 0);
        const auto itBegin = vPrefSum.begin();
        const auto itEnd = vPrefSum.begin() + idx;
        auto itSearch = std::lower_bound(itBegin,itEnd, diff);
        if (itSearch != itEnd)
        {
            if ((vPrefSum[idx] - *itSearch) > maxSum.maxSum)
            {
                maxSum.maxSum = vPrefSum[idx] - *itSearch;
                maxSum.startIdx = std::distance(itBegin, itSearch) + 1;
                maxSum.endIdx = idx;
            }
        }
        if (maxSum.maxSum == nCandies)
        {
            break;
        }
    }
}

void computeMaxSequenceParallel(const std::vector<int32_t>& vHomeCandies, const std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                                  MaxSum& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    bool found = false;
    std::vector<MaxSum> vMaxSums;

    uint32_t numThreads = omp_get_num_procs();
    omp_set_num_threads(numThreads);

    vMaxSums.resize(numThreads);
    #pragma  omp parallel for
    for (int32_t idx = 1; idx <= nHomes; idx ++)
    {
        if (!found)
        {
            uint32_t threadId = omp_get_thread_num();
            int32_t diff = std::max(vPrefSum[idx] - nCandies, 0);
            const auto itBegin = vPrefSum.begin();
            const auto itEnd = vPrefSum.begin() + idx;
            auto itSearch = std::lower_bound(itBegin,itEnd, diff);
            if (itSearch != itEnd)
            {
                if ((vPrefSum[idx] - *itSearch) > maxSum.maxSum)
                {
                    vMaxSums[threadId].maxSum = vPrefSum[idx] - *itSearch;
                    vMaxSums[threadId].startIdx = std::distance(itBegin, itSearch) + 1;
                    vMaxSums[threadId].endIdx = idx;
                }
            }
            if (vMaxSums[threadId].maxSum == nCandies)
            {
                found = true;
            }
        }
    }

    auto itSol = std::max_element(vMaxSums.begin(), vMaxSums.end());
    maxSum = *itSol;
}


int main() {
    int32_t nCandies;
    std::vector<int32_t> vHomeCandies;
    std::vector<int32_t> vPrefSum;
    MaxSum maxSum;;

    readFile(nCandies, vHomeCandies);
    computePrefixSums(vHomeCandies, vPrefSum);
    if (IMPL == 0)
    {
        computeMaxSequenceSequential(vHomeCandies, vPrefSum, nCandies, maxSum);
    }
    else if (IMPL == 1)
    {
        computeMaxSequenceParallel(vHomeCandies, vPrefSum, nCandies, maxSum);
    }
    else if (IMPL == 2)
    {

    }

    if (maxSum.maxSum != -1)
    {
        std::cout << "Start at home " << maxSum.startIdx << " and go to home "
        << maxSum.endIdx << " getting " << maxSum.maxSum << " pieces of candy" << std::endl;
    }
    else
    {
        std::cout << "Don't go here" << std::endl;
    }
    return 0;
}
