"""
Given a string s, find the length of the longest substring without repeating characters.



Example 1:

Input: s = "abcabcbb"
Output: 3
Explanation: The answer is "abc", with the length of 3.

Example 2:

Input: s = "bbbbb"
Output: 1
Explanation: The answer is "b", with the length of 1.

Example 3:

Input: s = "pwwkew"
Output: 3
Explanation: The answer is "wke", with the length of 3.
Notice that the answer must be a substring, "pwke" is a subsequence and not a substring.



Constraints:

    0 <= s.length <= 5 * 104
    s consists of English letters, digits, symbols and spaces.


"""
class Solution_me(object):
    def lengthOfLongestSubstring(self, s):
        """
        :type s: str
        :rtype: int
        """
        sub = []
        max = 0
        start_pos = 0
        i=0
        for c in s:
            if c not in sub:
                sub.append(c)
            else:
                sub.pop(0)
                if c not in sub:
                    sub.append(c)
            if len(sub) > max:
                max = len(sub)
                start_pos = i
            i+=1
        return max

test = Solution_me()
s = "pwwkew"
s = "bbbbb"
s = "abcabcbb"
s = "chicken"
s = "batty"
s = "foxie"
print(test.lengthOfLongestSubstring(s))
