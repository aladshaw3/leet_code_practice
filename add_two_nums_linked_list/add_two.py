# Definition for singly-linked list.
class ListNode(object):
    def __init__(self, val=0, next=None):
        self.val = val
        self.next = next

class Solution(object):
    def addTwoNumbers_chris(self, l1, l2):
        """
        :type l1: ListNode
        :type l2: ListNode
        :rtype: ListNode
        """
        first = 0
        i = l1
        j = 0
        while i is not None:
            first += pow(10, j) * i.val
            i = i.next
            j += 1

        second = 0
        i = l2
        j = 0
        while i is not None:
            second += pow(10, j) * i.val
            i = i.next
            j += 1

        result_integer = first + second

        j = 0
        result = None
        prev = None
        while result_integer != 0:
            tmp = ListNode(result_integer % 10)
            i = tmp
            if result is None:
                result = tmp
                prev = result

            prev.next = tmp
            prev = tmp
            result_integer = int(result_integer / 10)
            j += 1

        return result

    def addTwoNumbers_me_the_idiot(self, l1, l2):
        """
        :type l1: ListNode
        :type l2: ListNode
        :rtype: ListNode
        """
        res = []
        l1_list = []
        l1_list.append(l1.val)
        while (l1.next != None):
            l1 = l1.next
            l1_list.append(l1.val)

        l2_list = []
        l2_list.append(l2.val)
        while (l2.next != None):
            l2 = l2.next
            l2_list.append(l2.val)

        if len(l1_list) > len(l2_list):
            long_list = l1_list
            short_list = l2_list
        else:
            long_list = l2_list
            short_list = l1_list

        carry = 0
        i = 0
        res_list = []
        for v in long_list:
            sum = carry
            if i < len(short_list):
                sum += v + short_list[i]
            else:
                sum += v
            if sum >= 10:
                sum = 10 - sum
                carry = 1
            else:
                carry = 0
            res_list.append(sum)
            i+=1

        i = 0
        for r in res_list:
            res.append(ListNode(r,None))
            if i > 0 and i <= len(res_list):
                res[i-1].next = res[i]
            i+=1
        return res[0]

a = ListNode()
b = ListNode()
c = ListNode()

a.val = 1
a.next = b

b.val = 2
b.next = c

c.val = 3
c.next = None


d = ListNode()
e = ListNode()

d.val = 1
d.next = e

e.val = 8
e.next = None

# pass a and d
#   a actual num = 321
#   d actual num =  81
#       res =      402


res = Solution()
ll = res.addTwoNumbers_me_the_idiot(a,d)
print(ll.val)
while (ll.next != None):
    ll = ll.next
    print(ll.val)
