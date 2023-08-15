#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>

struct PrefSumIdx
{
    int32_t sum;
    int32_t val;
    PrefSumIdx(int32_t _sum, int32_t _val): sum(_sum), val(_val) {}

    bool operator<(const PrefSumIdx& other) const
    {
        return sum < other.sum;
    }
};

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

void computeMaxSeq(const std::vector<int32_t>& vHomeCandies, const std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                   int32_t& startIdx, int32_t& endIdx, int32_t& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    std::set<PrefSumIdx> sPrefSums;

    maxSum = -1;
    startIdx = -1;
    endIdx = -1;

    sPrefSums.insert({vPrefSum[0], 0});
    for (int32_t idx = 1; idx <= nHomes; idx ++)
    {
        int32_t diff = std::max(vPrefSum[idx] - nCandies, 0);
        auto itSearch = sPrefSums.lower_bound({diff, 0});
        if (itSearch != sPrefSums.end())
        {
            if ((vPrefSum[idx] - itSearch->sum) > maxSum)
            {
                maxSum = vPrefSum[idx] - itSearch->sum;
                startIdx = itSearch->val + 1;
                endIdx = idx;
            }
        }
        if (maxSum == nCandies)
        {
            break;
        }
        sPrefSums.insert({vPrefSum[idx],idx});
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
    computeMaxSeq(vHomeCandies, vPrefSum, nCandies, startIdx, endIdx, maxSum);

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
