"""
Given two sorted arrays nums1 and nums2 of size m and n respectively, return the median of the two sorted arrays.

The overall run time complexity should be O(log (m+n)).



Example 1:

Input: nums1 = [1,3], nums2 = [2]
Output: 2.00000
Explanation: merged array = [1,2,3] and median is 2.

Example 2:

Input: nums1 = [1,2], nums2 = [3,4]
Output: 2.50000
Explanation: merged array = [1,2,3,4] and median is (2 + 3) / 2 = 2.5.



Constraints:

    nums1.length == m
    nums2.length == n
    0 <= m <= 1000
    0 <= n <= 1000
    1 <= m + n <= 2000
    -106 <= nums1[i], nums2[i] <= 106

"""

class Solution(object):
    def findMedianSortedArrays(self, nums1, nums2):
        """
        :type nums1: List[int]
        :type nums2: List[int]
        :rtype: float
        """
        if len(nums1) > len(nums2):
            long_list = nums1
            short_list = nums2
        else:
            long_list = nums2
            short_list = nums1
        total_length = len(long_list) + len(short_list)
        rem = total_length % 2
        double_index = False
        if rem == 0:
            double_index = True
        start_index = int((total_length-1)/2)

        if double_index == True:
            r = start_index+1
        else:
            r = start_index

        merged = []
        llind=0
        slind=0
        for i in range(r+1):
            if slind < len(short_list):
                if long_list[llind] < short_list[slind]:
                    merged.append(long_list[llind])
                    llind+=1
                else:
                    merged.append(short_list[slind])
                    slind+=1
            else:
                merged.append(long_list[llind])
                llind+=1

        if double_index == False:
            return merged[start_index]
        else:
            return (merged[start_index]+merged[start_index+1])/2

    def findMedianSortedArrays_BinarySearch(self, nums1, nums2):
        """
        :type nums1: List[int]
        :type nums2: List[int]
        :rtype: float
        """
        A, B = nums1, nums2
        total = len(nums1) + len(nums2)
        half = total // 2  # use // for int divide

        if len(B) < len(A):
            A, B = B, A

        # log( min(n,m) )
        l, r = 0, len(A) - 1
        while True:
            i = (l+r) // 2 #A
            j = half - i - 2 #B
            Aleft = A[i] if i >=0 else float("-infinity")
            Aright = A[i+1] if (i+1) < len(A) else float("infinity")
            Bleft = B[j] if j >= 0 else float("-infinity")
            Bright = B[j+1] if (j+1) < len(B) else float("infinity")

            # partition is correct
            if Aleft <= Bright and Bleft <= Aright:
                # odd
                if total % 2:
                    return min(Aright, Bright)
                # even
                return (max(Aleft,Bleft) + min(Aright,Bright)) / 2
            elif Aleft > Bright:
                r = i - 1
            else:
                l = i + 1


nums1 = [1,2]
nums2 = [3]
test = Solution()
print(test.findMedianSortedArrays_BinarySearch(nums1,nums2))
