#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <omp.h>
#include <cassert>
#include <chrono>

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

    // Required for the reduction step in the parallelization.
    bool operator<(const MaxSum& other) const
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

/**
 * A function that reads specification from the file named fileInputFrom.
 * @param nCandies returns the maximal number of candies a user can take.
 * @param vHomeCandies returns the vector of candies each home gives in sequential order
 * @param fileInputName passes the file name.
 */
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
/**
 * The computation of a prefix sums of number of candies each home gives starting from 0 up to that home.
 * We assume home 0 gives 0 for the sake of brevity of the code.
 * @param vHomeCandies the vector of candies each home gives in sequential order
 * @param vPrefSum returns the vector that contains the prefix sum.
 */
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

/**
 * Parallelized computation prefix sums of number of candies each home gives starting from 0 up to that home.
 * The method computes by exploiting the rule that
 * prefSum(start, end) = prefSum(start, mid) + prefSum(mid+1, end) where prefSum(i, j) denotes a sum of the candies
 * starting at the home i and sequentially all the way up to end.
 * 1) First we divide homes in chunks and compute prefSum(start, end) for each of the chunks by using classical
 * method to compute prefixSum.
 * 2) Then we add these computed sums to the next chunk incrementally.
 * Time complexity is O(n / k) where k is a number rof processors, although the method is limited by the memory
 * bandwidth.
 * We assume home 0 gives 0 for the sake of brevity of the code.
 * @param vHomeCandies the vector of candies each home gives in sequential order
 * @param vPrefSum returns the vector that contains the prefix sum.
 */
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

    #pragma  omp parallel  for schedule(static) default(none) shared(nHomes, vSums, vHomeCandies, vPrefSum)
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

    #pragma  omp parallel  for schedule(static) default(none) shared(nHomes, vSums, vHomeCandies, vPrefSum)
    for (uint32_t idx = 0; idx < nHomes; idx++)
    {
        uint32_t threadId = omp_get_thread_num();
        if (threadId > 0)
        {
            vPrefSum[idx+1] += vSums[threadId-1];
        }
    }
}

/**
 * A classical sequential search where for each home we check what would be the starting home
 * that maximizes the number of candies one child can get. We use binary search for this
 * to reduce the time complexity of the algorithm from O(n^2) to O(nlogn).
 *
 * @param vHomeCandies the vector of candies each home gives in sequential order
 * @param vPrefSum the vector that contains the prefix sum.
 * @param nCandies the maximum number of candies each child can get
 * @param maxSum returns a maximum number of candies that a child can get alongside the starting idx and end idx.
 */
void computeMaxSequenceSequential(const std::vector<int32_t>& vHomeCandies, const std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                   MaxSum& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    maxSum = MaxSum();
    auto itBegin = vPrefSum.begin();

    for (int32_t idx = 1; idx <= nHomes; idx ++)
    {
        int32_t diff = std::max(vPrefSum[idx] - nCandies, 0);

        const auto itEnd = vPrefSum.begin() + idx;
        auto itSearch = std::lower_bound(itBegin,itEnd, diff);
        if (itSearch != itEnd)
        {
            int32_t curMaxSum = vPrefSum[idx] - *itSearch;
            // If the curMaxSum == maxSum.maxSum, the start Idx will be for greater or equal to
            // the startIdx of the currently found maxSum since all elements are positive.
            // thus we should not consider this case.
            if (curMaxSum > maxSum.maxSum)
            {
                maxSum.maxSum = curMaxSum;
                maxSum.startIdx = std::distance(vPrefSum.begin(), itSearch) + 1;
                maxSum.endIdx = idx;
            }
            itBegin = itSearch;
        }

        if (maxSum.maxSum == nCandies)
        {
            break;
        }
    }
}

/**
 * A classical parallelized sequential search where for each home we check what would be the starting home
 * that maximizes the number of candies one child can get. We use binary search for this
 * to reduce the time complexity of the algorithm from O(n^2 / k) to O(nlogn / k), where k is
 * a number of threads. For a parallelization, we use map reduce where in the reduction step we use stl algorithm
 * to fin maximum of already computed sublocal maximums.
 *
 * @param vHomeCandies the vector of candies each home gives in sequential order
 * @param vPrefSum the vector that contains the prefix sum.
 * @param nCandies the maximum number of candies each child can get
 * @param maxSum returns a maximum number of candies that a child can get alongside the starting idx and end idx.
 */
void computeMaxSequenceParallel(const std::vector<int32_t>& vHomeCandies, std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                                  MaxSum& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    bool found = false;
    std::vector<MaxSum> vMaxSums;
    std::vector<std::vector<int32_t >::iterator > vBeginIterators;

    uint32_t numThreads = omp_get_num_procs();
    omp_set_num_threads(numThreads);
    vMaxSums.resize(numThreads);
    vBeginIterators.resize(numThreads);
    for (auto& it: vBeginIterators)
    {
        it = vPrefSum.begin();
    }

    #pragma  omp parallel  for schedule(static) default(none) shared(nHomes, found, vPrefSum, nCandies, vMaxSums, vBeginIterators)
    for (int32_t idx = 1; idx <= nHomes; idx ++)
    {
        if (!found)
        {
            uint32_t threadId = omp_get_thread_num();
            int32_t diff = std::max(vPrefSum[idx] - nCandies, 0);
            const auto itBegin = vBeginIterators[threadId];
            const auto itEnd = vPrefSum.begin() + idx;
            auto itSearch = std::lower_bound(itBegin,itEnd, diff);
            if (itSearch != itEnd)
            {
                int32_t curMaxSum = vPrefSum[idx] - *itSearch;
                if (curMaxSum > vMaxSums[threadId].maxSum)
                {
                    vMaxSums[threadId].maxSum = curMaxSum;
                    vMaxSums[threadId].startIdx = std::distance(vPrefSum.begin(), itSearch) + 1;
                    vMaxSums[threadId].endIdx = idx;
                }
                vBeginIterators[threadId] = itSearch;
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


/**
 * A window search where we keep bounds of window which currently contains  the starting home
 * that maximizes the number of candies one child can get for the currently ending chome.
 * The time complexity of this approach is O(n) since each iterator can move at most n times.
 *
 * @param vHomeCandies the vector of candies each home gives in sequential order
 * @param vPrefSum the vector that contains the prefix sum.
 * @param nCandies the maximum number of candies each child can get
 * @param maxSum returns a maximum number of candies that a child can get alongside the starting idx and end idx.
 */
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


/**
* A parallelized window search where we keep bounds of window which currently contains  the starting home
        * that maximizes the number of candies one child can get for the currently ending chome.
* The time complexity of this approach is O(n) since the beginning iterator can iterate O(n) times, but
 * in general case it can bring improvements up to O(n / k) where k is a number of threads.
*
* @param vHomeCandies the vector of candies each home gives in sequential order
* @param vPrefSum the vector that contains the prefix sum.
* @param nCandies the maximum number of candies each child can get
* @param maxSum returns a maximum number of candies that a child can get alongside the starting idx and end idx.
*/
void computeMaxSequenceWindowApproachParallel(const std::vector<int32_t>& vHomeCandies, const std::vector<int32_t>& vPrefSum, int32_t& nCandies,
                                      MaxSum& maxSum)
{
    const uint32_t nHomes = vHomeCandies.size();
    std::vector<MaxSum> vMaxSums;
    std::vector<bool> vInit;

    std::vector<int32_t > vStartIdx;
    uint32_t numThreads = omp_get_num_procs();
    omp_set_num_threads(numThreads);
    vStartIdx.resize(numThreads);
    vMaxSums.resize(numThreads);
    vInit.resize(numThreads);

    #pragma  omp parallel  for schedule(static) default(none) shared(nHomes, vPrefSum, nCandies, vMaxSums, vInit, vStartIdx)
    for (int32_t endIdx = 1; endIdx <= nHomes; endIdx ++)
    {
        uint32_t threadId = omp_get_thread_num();
        if (!vInit[threadId])
        {
            int32_t diff = std::max(vPrefSum[endIdx] - nCandies, 0);;
            const auto itEnd = vPrefSum.begin() + endIdx;
            auto itSearch = std::lower_bound(vPrefSum.begin(),itEnd, diff);
            vStartIdx[threadId] = std::distance(vPrefSum.begin(), itSearch);
            vInit[threadId] = true;
        }
        while (vPrefSum[endIdx]  - vPrefSum[vStartIdx[threadId]] > nCandies)
        {
            vStartIdx[threadId]++;
        }
        if (endIdx == vStartIdx[threadId])
        {
            continue;
        }
        int32_t curMaxSum = vPrefSum[endIdx] - vPrefSum[vStartIdx[threadId]];
        if (curMaxSum > vMaxSums[threadId].maxSum)
        {
            vMaxSums[threadId].maxSum = curMaxSum;
            vMaxSums[threadId].startIdx = vStartIdx[threadId]+1;
            vMaxSums[threadId].endIdx = endIdx;
        }
    }
    auto itSol = std::max_element(vMaxSums.begin(), vMaxSums.end());
    maxSum = *itSol;
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
    else if (implAlg == 3)
    {
        computeMaxSequenceWindowApproachParallel(vHomeCandies, vPrefSum, nCandies, maxSum);
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

    checkHomes(fileInputName, maxSum, 3, 1);
    EXPECT_EQ(maxSum.maxSum, 10);
    EXPECT_EQ(maxSum.startIdx, 2);
    EXPECT_EQ(maxSum.endIdx, 5);
}


void notEnterTest()
{
    const std::string fileInputName = "tests/input_not_enter.txt";
    MaxSum maxSum;

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

    checkHomes(fileInputName, maxSum, 3, 1);
    EXPECT_EQ(maxSum.maxSum, -1);
    EXPECT_EQ(maxSum.startIdx, -1);
    EXPECT_EQ(maxSum.endIdx, -1);

}

void altTest(void)
{
    const std::string fileInputName = "tests/input_alt.txt";
    MaxSum maxSum;

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

    checkHomes(fileInputName, maxSum, 3, 1);
    EXPECT_EQ(maxSum.maxSum, 9);
    EXPECT_EQ(maxSum.startIdx, 8);
    EXPECT_EQ(maxSum.endIdx, 10);
}

void zeroTest(void)
{
    const std::string fileInputName = "tests/input_zero.txt";
    MaxSum maxSum;

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

    checkHomes(fileInputName, maxSum, 3, 1);
    EXPECT_EQ(maxSum.maxSum, 0);
    EXPECT_EQ(maxSum.startIdx, 5);
    EXPECT_EQ(maxSum.endIdx, 5);
}

void zeroNotTest(void)
{
    const std::string fileInputName = "tests/input_zero_not.txt";
    MaxSum maxSum;

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

    checkHomes(fileInputName, maxSum, 3, 1);
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
    //MaxSum maxSum;
    //checkHomes("input.txt", maxSum, 1, 1);
    return 0;
}
