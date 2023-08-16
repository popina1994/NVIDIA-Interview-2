#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <omp.h>
#include <cassert>

#define EXPECT_EQ(A, B) assert(A == B);

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
        // Also required due to the condition specified in the problem.
        if (maxSum == other.maxSum)
        {
            return startIdx > other.startIdx;
        }
        return maxSum < other.maxSum;
    }
};

void readFile(int32_t& nCandies, std::vector<int32_t>& vHomeCandies, const std::string& fileInputName)
{
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

void computePrefixSumsParallel(const std::vector<int32_t>& vHomeCandies, std::vector<int32_t>& vPrefSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    std::vector<int32_t> vSums;
    vPrefSum.resize(nHomes+1);
    int32_t curSum = 0;
    vPrefSum[0] = 0;

    uint32_t numThreads = omp_get_num_procs();
    omp_set_num_threads(numThreads);
    vSums.resize(numThreads);

    #pragma  omp parallel for default(none) shared(nHomes, vSums, vHomeCandies, vPrefSum)
    for (uint32_t idx = 0; idx < nHomes; idx++)
    {
        uint32_t threadId = omp_get_thread_num();
        vSums[threadId] += vHomeCandies[idx];
        vPrefSum[idx+1] = vSums[threadId];
    }

    for (uint32_t idx = 1; idx < numThreads; idx++)
    {
        vSums[idx] += vSums[idx-1];
    }

    #pragma  omp parallel for default(none) shared(nHomes, vSums, vHomeCandies, vPrefSum)
    for (uint32_t idx = 0; idx < nHomes; idx++)
    {
        uint32_t threadId = omp_get_thread_num();
        if (threadId > 0)
        {
            vPrefSum[idx+1] += vSums[threadId-1];
        }
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
            int32_t curMaxSum = vPrefSum[idx] - *itSearch;
            if (curMaxSum > maxSum.maxSum)
            {
                maxSum.maxSum = curMaxSum;
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

    #pragma  omp parallel for default(none) shared(nHomes, found, vPrefSum, nCandies, vMaxSums)
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
                int32_t curMaxSum = vPrefSum[idx] - *itSearch;
                if (curMaxSum > vMaxSums[threadId].maxSum)
                {
                    vMaxSums[threadId].maxSum = curMaxSum;
                    vMaxSums[threadId].startIdx = std::distance(itBegin, itSearch) + 1;
                    vMaxSums[threadId].endIdx = idx;
                }
            }
            if ((threadId == 0) && (vMaxSums[threadId].maxSum == nCandies))
            {
                found = true;
            }
        }
    }

    auto itSol = std::max_element(vMaxSums.begin(), vMaxSums.end());
    maxSum = *itSol;
}


void computeMaxSequenceWindowApproach(const std::vector<int32_t>& vHomeCandies, const std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                                  MaxSum& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    maxSum = MaxSum();
    int32_t startIdx = 0;

    for (int32_t endIdx = 1; endIdx <= nHomes; endIdx ++)
    {
        while (vPrefSum[endIdx]  - vPrefSum[startIdx] > nCandies)
        {
            startIdx++;
        }
        if (endIdx == startIdx)
        {
            continue;
        }
        int32_t curMaxSum = vPrefSum[endIdx] - vPrefSum[startIdx];
        if (curMaxSum > maxSum.maxSum)
        {
            maxSum.maxSum = curMaxSum;
            maxSum.startIdx = startIdx+1;
            maxSum.endIdx = endIdx;
        }
        if (maxSum.maxSum == nCandies)
        {
            break;
        }
    }
}


void checkHomes(const std::string& fileInputName, MaxSum& maxSum, int32_t implAlg, int32_t prefSumImpl)
{
    int32_t nCandies;
    std::vector<int32_t> vHomeCandies;
    std::vector<int32_t> vPrefSum;

    readFile(nCandies, vHomeCandies, fileInputName);
    if (prefSumImpl == 0)
    {
        computePrefixSums(vHomeCandies, vPrefSum);
    }
    else
    {
        computePrefixSumsParallel(vHomeCandies, vPrefSum);
    }

    if (implAlg == 0)
    {
        computeMaxSequenceSequential(vHomeCandies, vPrefSum, nCandies, maxSum);
    }
    else if (implAlg == 1)
    {
        computeMaxSequenceParallel(vHomeCandies, vPrefSum, nCandies, maxSum);
    }
    else if (implAlg == 2)
    {
        computeMaxSequenceWindowApproach(vHomeCandies, vPrefSum, nCandies, maxSum);
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
}

void defaultTest()
{
    const std::string fileInputName = "tests/input_default.txt";
    MaxSum maxSum;
    maxSum = {-1, -1, -1};
    checkHomes(fileInputName, maxSum, 0, 1);
    EXPECT_EQ(maxSum.maxSum, 10);
    EXPECT_EQ(maxSum.startIdx, 2);
    EXPECT_EQ(maxSum.endIdx, 5);

    checkHomes(fileInputName, maxSum, 1, 1);
    EXPECT_EQ(maxSum.maxSum, 10);
    EXPECT_EQ(maxSum.startIdx, 2);
    EXPECT_EQ(maxSum.endIdx, 5);

    checkHomes(fileInputName, maxSum, 2, 1);
    EXPECT_EQ(maxSum.maxSum, 10);
    EXPECT_EQ(maxSum.startIdx, 2);
    EXPECT_EQ(maxSum.endIdx, 5);
}


void notEnterTest()
{
    const std::string fileInputName = "tests/input_not_enter.txt";
    MaxSum maxSum;
    maxSum = {-1, -1, -1};
    checkHomes(fileInputName, maxSum, 0, 0);
    EXPECT_EQ(maxSum.maxSum, -1);
    EXPECT_EQ(maxSum.startIdx, -1);
    EXPECT_EQ(maxSum.endIdx, -1);

    checkHomes(fileInputName, maxSum, 1, 0);
    EXPECT_EQ(maxSum.maxSum, -1);
    EXPECT_EQ(maxSum.startIdx, -1);
    EXPECT_EQ(maxSum.endIdx, -1);

    checkHomes(fileInputName, maxSum, 2, 0);
    EXPECT_EQ(maxSum.maxSum, -1);
    EXPECT_EQ(maxSum.startIdx, -1);
    EXPECT_EQ(maxSum.endIdx, -1);

}

void altTest(void)
{
    const std::string fileInputName = "tests/input_alt.txt";
    MaxSum maxSum;
    maxSum = {-1, -1, -1};
    checkHomes(fileInputName, maxSum, 0, 1);
    EXPECT_EQ(maxSum.maxSum, 9);
    EXPECT_EQ(maxSum.startIdx, 8);
    EXPECT_EQ(maxSum.endIdx, 10);

    checkHomes(fileInputName, maxSum, 1, 1);
    EXPECT_EQ(maxSum.maxSum, 9);
    EXPECT_EQ(maxSum.startIdx, 8);
    EXPECT_EQ(maxSum.endIdx, 10);

    checkHomes(fileInputName, maxSum, 2, 1);
    EXPECT_EQ(maxSum.maxSum, 9);
    EXPECT_EQ(maxSum.startIdx, 8);
    EXPECT_EQ(maxSum.endIdx, 10);
}

void zeroTest(void)
{
    const std::string fileInputName = "tests/input_zero.txt";
    MaxSum maxSum;
    maxSum = {-1, -1, -1};
    checkHomes(fileInputName, maxSum, 0, 1);
    EXPECT_EQ(maxSum.maxSum, 0);
    EXPECT_EQ(maxSum.startIdx, 5);
    EXPECT_EQ(maxSum.endIdx, 5);

    checkHomes(fileInputName, maxSum, 1, 1);
    EXPECT_EQ(maxSum.maxSum, 0);
    EXPECT_EQ(maxSum.startIdx, 5);
    EXPECT_EQ(maxSum.endIdx, 5);

    checkHomes(fileInputName, maxSum, 2, 1);
    EXPECT_EQ(maxSum.maxSum, 0);
    EXPECT_EQ(maxSum.startIdx, 5);
    EXPECT_EQ(maxSum.endIdx, 5);
}

void zeroNotTest(void)
{
    const std::string fileInputName = "tests/input_zero_not.txt";
    MaxSum maxSum;
    maxSum = {-1, -1, -1};
    checkHomes(fileInputName, maxSum, 0, 1);
    EXPECT_EQ(maxSum.maxSum, -1);
    EXPECT_EQ(maxSum.startIdx, -1);
    EXPECT_EQ(maxSum.endIdx, -1);

    checkHomes(fileInputName, maxSum, 1, 1);
    EXPECT_EQ(maxSum.maxSum, -1);
    EXPECT_EQ(maxSum.startIdx, -1);
    EXPECT_EQ(maxSum.endIdx, -1);

    checkHomes(fileInputName, maxSum, 2, 1);
    EXPECT_EQ(maxSum.maxSum, -1);
    EXPECT_EQ(maxSum.startIdx, -1);
    EXPECT_EQ(maxSum.endIdx, -1);
}


void runTests()
{
    defaultTest();
    notEnterTest();
    altTest();
    zeroTest();
    zeroNotTest();
}

int main()
{
    runTests();
    return 0;
}
