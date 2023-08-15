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
    vPrefSum.resize(nHomes);
    int32_t curSum = 0;
    for (uint32_t idx = 0; idx < nHomes; idx++)
    {
        curSum += vHomeCandies[idx];
        vPrefSum[idx] = curSum;
    }
}

int main() {
    int32_t nCandies;
    std::vector<int32_t> vHomeCandies;
    std::vector<int32_t> vPrefSum;

    readFile(nCandies, vHomeCandies);
    computePrefixSums(vHomeCandies, vPrefSum);

    const uint32_t nHomes = vHomeCandies.size();
    std::set<PrefSumIdx> sPrefSums;

    int32_t maxSum = -1;
    int32_t startIdx = -1;
    int32_t endIdx = -1;

    for (int32_t idx = 0; idx < nHomes; idx ++)
    {
        if (nCandies - vHomeCandies[idx] >= 0)
        {
            maxSum = vHomeCandies[idx];
            startIdx = idx;
            endIdx = idx;
        }
        auto itSearch = sPrefSums.lower_bound({nCandies - vPrefSum[idx], 0});
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

    if (maxSum != -1)
    {
        std::cout << "Start at home " << startIdx + 1 << " and go to home "
        << endIdx + 1 << " getting " << maxSum << " pieces of candy" << std::endl;
    }
    else
    {
        std::cout << "Don't go here" << std::endl;
    }
    return 0;
}
