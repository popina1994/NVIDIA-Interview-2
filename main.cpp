#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>

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
                   int32_t& startIdx, int32_t& endIdx, int32_t& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();

    maxSum = -1;
    startIdx = -1;
    endIdx = -1;

    for (int32_t idx = 1; idx <= nHomes; idx ++)
    {
        int32_t diff = std::max(vPrefSum[idx] - nCandies, 0);
        auto itSearch = std::lower_bound(vPrefSum.begin(),
                                            vPrefSum.begin() + idx, diff);
        if (itSearch != vPrefSum.end())
        {
            if ((vPrefSum[idx] - *itSearch) > maxSum)
            {
                maxSum = vPrefSum[idx] - *itSearch;
                startIdx = std::distance(vPrefSum.begin(), itSearch) + 1;
                endIdx = idx;
            }
        }
        if (maxSum == nCandies)
        {
            break;
        }
    }
}

void computeMaxSequenceParallel(const std::vector<int32_t>& vHomeCandies, const std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                                  int32_t& startIdx, int32_t& endIdx, int32_t& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();

    maxSum = -1;
    startIdx = -1;
    endIdx = -1;

    #pragma  omp parallel for
    for (int32_t idx = 1; idx <= nHomes; idx ++)
    {
        int32_t diff = std::max(vPrefSum[idx] - nCandies, 0);
        auto itSearch = std::lower_bound(vPrefSum.begin(),
                                         vPrefSum.begin() + idx, diff);
        if (itSearch != vPrefSum.end())
        {
            if ((vPrefSum[idx] - *itSearch) > maxSum)
            {
                maxSum = vPrefSum[idx] - *itSearch;
                startIdx = std::distance(vPrefSum.begin(), itSearch) + 1;
                endIdx = idx;
            }
        }
        if (maxSum == nCandies)
        {
            break;
        }
    }
}


int main() {
    int32_t nCandies;
    std::vector<int32_t> vHomeCandies;
    std::vector<int32_t> vPrefSum;
    int32_t maxSum;
    int32_t startIdx;
    int32_t endIdx;

    readFile(nCandies, vHomeCandies);
    computePrefixSums(vHomeCandies, vPrefSum);
    computeMaxSequenceSequential(vHomeCandies, vPrefSum, nCandies, startIdx, endIdx, maxSum);

    if (maxSum != -1)
    {
        std::cout << "Start at home " << startIdx << " and go to home "
        << endIdx << " getting " << maxSum << " pieces of candy" << std::endl;
    }
    else
    {
        std::cout << "Don't go here" << std::endl;
    }
    return 0;
}
