'''
Given an array of integers nums and an integer target, return indices of the two numbers such that they add up to target.

You may assume that each input would have exactly one solution, and you may not use the same element twice.

You can return the answer in any order.

Example 1:

Input: nums = [2,7,11,15], target = 9
Output: [0,1]
Explanation: Because nums[0] + nums[1] == 9, we return [0, 1].

Example 2:

Input: nums = [3,2,4], target = 6
Output: [1,2]

Example 3:

Input: nums = [3,3], target = 6
Output: [0,1]



Constraints:

    2 <= nums.length <= 10^4
    -10^9 <= nums[i] <= 10^9
    -10^9 <= target <= 10^9
    Only one valid answer exists.

'''

class Solution(object):
    def twoSum_via_hash(self, nums, target):
        """
        :type nums: List[int]
        :type target: int
        :rtype: List[int]

        Time complexity O(n)
        Space complexity O(n)

        How it works...
            Dynamically build a hash table of prior values
            observed. Start by calculating the value needed
            based on the current iterate value 'n' and the target.
            Search the hashtable for the needed value 'm'. If not
            present, store current iterate value in hashtable and
            continue.

            The hashtable grows with each iteration and will be at
            most O(n) in size.
        """
        hashmap = {}
        for i, n in enumerate(nums):
            m = target - n
            if m in hashmap:
                return [hashmap[m], i]
            hashmap[n] = i

    def twoSum_via_bruteforce(self, nums, target):
        """
        :type nums: List[int]
        :type target: int
        :rtype: List[int]

        Time complexity O(n^2)
        Space complexity O(1)

	How it works...
	    Iterate through the loop 2x times in an outer and
            inner loop. Check the sum and report.
        """
        for i in range(len(nums)):
            for j in range(i + 1, len(nums)):
                if nums[j] == target - nums[i]:
                    return [i, j]

obj = Solution()
print(obj.twoSum_via_hash([1,2,3,4,5,6],7))
print(obj.twoSum_via_bruteforce([1,2,3,4,5,6],7))
