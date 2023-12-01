from Vector import *

# Program
DP1 = {}
def edit_distance(rna1, rna2):
    """ This edit-distance function takes the input of two strings. Its output is a dynamic-programming table, where each indices is a tuple that includes the following information: current edit-distance,
    the count of different approaches to achieve edit-distance, the operations of that index. 
    """
    len1 = len(rna1)
    len2 = len(rna2)
    
    # set up three empty tables to seperately record DP table for edit-distance, DP-count table for number of different approaches, DP_ops table for all the corresponding operations
    DP_ = [[0 for i in range(len2 + 1)] for j in range(len1 + 1)]
    DP_count = [[0 for i in range(len2 + 1)] for j in range(len1 + 1)]
    DP_ops = [[[] for _ in range(len2 + 1)] for _ in range(len1 + 1)]   
    
    # calculate the table using dynamic programming 
    for i in range(len1 + 1):
        for j in range(len2 + 1):
            if i == 0 and j == 0: 
                DP_[i][j] = j
                DP_count[i][j] = 1
            elif i == 0:
                DP_[i][j] = j
                DP_count[i][j] = 1
                DP_ops[i][j].append((0, -1))
            elif j == 0:
                DP_[i][j] = i
                DP_count[i][j] = 1
                DP_ops[i][j].append((-1, 0))
            elif rna1[i - 1] == rna2[j - 1]:
                DP_[i][j] = DP_[i - 1][j - 1]
                DP_count[i][j] = DP_count[i - 1][j - 1]
                DP_ops[i][j].append((-1, -1))
            else:
                DP_[i][j] = 1 + min(DP_[i - 1][j], DP_[i][j - 1], DP_[i - 1][j - 1])
                
                # Count all the ways to reach the min edit distance
                # Update the operations table 
                if DP_[i][j] == 1 + DP_[i - 1][j]:
                    DP_count[i][j] += DP_count[i - 1][j]
                    DP_ops[i][j].append((-1, 0))
                if DP_[i][j] == 1 + DP_[i][j - 1]:
                    DP_count[i][j] += DP_count[i][j - 1]
                    DP_ops[i][j].append((0, -1))
                if DP_[i][j] == 1 + DP_[i - 1][j - 1]:
                    DP_count[i][j] += DP_count[i - 1][j - 1]
                    DP_ops[i][j].append((-1, -1))
            DP1.update({(i,j): [DP_[i][j], DP_count[i][j], DP_ops[i][j]]})
    return DP1


def diam(a,b,a_,b_, DP, memo = {}):
    """ This diameter function takes six inputs, four integers, one DP1 table dictionary we produced in the last function, and one empty DP2 memo to help us quickly sort the number. Its output is a single number, which is the diameter, maximum difference between 
    two paths. One path is from (a, b) to (0, 0) and the other path is from (a_, b_) to (0, 0).
    """

    # Check whether the memo is filled already
    if (a, b, a_, b_) in memo:
        return memo[(a, b, a_, b_)]
    
    result = 0
    options = []
    if a == b == a_ == b_ == 0:
        return 0
    # the special case of (a,b) and (a_,b_) are the same point
    if a == a_ and b == b_: 
        # if there is only one choice of the operation 
        if len(DP[a,b][2]) == 1:  
            result = diam(a + DP[a,b][2][0][0], b + DP[a,b][2][0][1], a_ + DP[a,b][2][0][0], b_ + DP[a,b][2][0][1], DP) 
            memo[(a, b, a_, b_)] = result
            return result 
        # if there is only two choices of the operation
        elif len(DP[a,b][2]) == 2:
            for i in range(2):
                for j in range(i,2):
                        (a1,b1) = DP[a,b][2][i]
                        (a2,b2) = DP[a_,b_][2][j]
                        option = diam(a + a1, b + b1, a_ + a2, b_ + b2, DP) + 2
                        options.append(option)
                        result = max(options)
            memo[(a, b, a_, b_)] = result
            return result 
        # if there is only three choices of the operation
        elif len(DP[a,b][2]) == 3:
            for i in range(3):
                for j in range(i,3):
                        (a1,b1) = DP[a,b][2][i]
                        (a2,b2) = DP[a_,b_][2][j]
                        if DP[a,b][2][i] == DP[a_,b_][2][j]:
                            option = diam(a + a1, b + b1, a_ + a2, b_ + b2, DP)
                            options.append(option)
                        else: 
                            option = diam(a + a1, b + b1, a_ + a2, b_ + b2, DP) + 2
                            options.append(option)
                        result = max(options)
            memo[(a, b, a_, b_)] = result
            return result 
    # (a,b) subsumes (a_,b_) and neither of them subsumes each other       
    elif (a > a_ and b > b_) or (a == a_ and b > b_) or (a > a_ and b == b_) or (a > a_ and b < b_) or (a < a_ and b > b_) :      
        for arrow in DP[a,b][2]:
            option = diam(a + arrow[0], b + arrow[1], a_, b_, DP) + 1
            result = max(result, option)
        memo[(a, b, a_, b_)] = result
        return result
    # (a_,b_) subsumes (a,b)
    else: 
        for arrow in DP[a_,b_][2]:
            option = diam(a, b, a_ + arrow[0], b_ + arrow[1], DP) + 1
            result = max(result, option)
        memo[(a, b, a_, b_)] = result
        return result



def dist_vector(x,y,x_,y_, DP, memo1 = {}):
    """This pairwise distance vector function aims to produce the pairwise distance for all possible pairs. It takes the input of four integers, one DP table dictionary we produced at the beginning, 
    and one empty memo1 to help us quickly sort the number. Its output is a list including a series of numbers, each stands for the number of pairs of paths that have such distance. The pair of path would be such that one path is from (x, y) to (0, 0) and the other path is from (x_, y_) to (0, 0).
    """
    # Initiate the solution by applying vector class
    solution = Vector(2*(len(rna1)+len(rna2)))

    if (x, y, x_, y_) in memo1:
        return memo1[(x, y, x_, y_)]
    
    if x == y == x_ == y_ == 0:
        solution.set(0,1)
        return solution
    # the special case of (x,y) and (x_,y_) are the same point
    if x == x_ and y == y_:  
        if len(DP[x,y][2]) == 1:
            solution = dist_vector(x + DP[x,y][2][0][0], y + DP[x,y][2][0][1], x_ + DP[x,y][2][0][0], y_ + DP[x,y][2][0][1], DP)
            memo1[(x, y, x_, y_)] = solution
            return solution 
            
        elif len(DP[x,y][2]) == 2:
            for i in range(2):
                for j in range(i,2):
                    (x1,y1) = DP[x,y][2][i]
                    (x2,y2) = DP[x_,y_][2][j]
                    if DP[x,y][2][i] == DP[x_,y_][2][j]:
                        print((x,y,x_,y_), dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP))
                        solution = solution + dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP).shift(0)
                    else: 
                        print((x,y,x_,y_), dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP))
                        solution = solution + dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP).shift(2)
            memo1[(x, y, x_, y_)] = solution
            return solution 
            
        elif len(DP[x,y][2]) == 3:
            for i in range(3):
                for j in range(i,3):
                    (x1,y1) = DP[x,y][2][i]
                    (x2,y2) = DP[x_,y_][2][j]
                    if DP[x,y][2][i] == DP[x_,y_][2][j]:
                        print((x + x1,y + y1,x_ + x2,y_ + y2), dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP))
                        solution = solution + dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP).shift(0)
                    else: 
                        print((x + x1,y + y1,x_ + x2,y_ + y2), dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP))
                        solution = solution + dist_vector(x + x1, y + y1, x_ + x2, y_ + y2, DP).shift(2)
            memo1[(x, y, x_, y_)] = solution
            return solution
    # (x,y) subsumes (x_,y_) and neither of them subsumes each other 
    elif (x > x_ and y > y_) or (x == x_ and y > y_) or (x > x_ and y == y_) or (x > x_ and y < y_) or (x < x_ and y > y_) :     
        for arrow in DP[x,y][2]:
            solution = solution + dist_vector(x + arrow[0], y + arrow[1], x_, y_, DP).shift(1)
        memo1[(x, y, x_, y_)] = solution
        return solution
            
    # (x_,y_) subsumes (x,y)
    else: 
        for arrow in DP[x_,y_][2]: 
            solution = solution + dist_vector(x, y, x_ + arrow[0], y_ + arrow[1], DP).shift(1)
        memo1[(x, y, x_, y_)] = solution
        return solution




rna1 = "AATT"
rna2 = "TTAA"
DP1 = edit_distance(rna1, rna2)
print(DP1)
print("diameter: ", diam(len(rna1),len(rna2),len(rna1),len(rna2),DP1))
result = dist_vector(len(rna1),len(rna2),len(rna1),len(rna2),DP1)
print("pairwise distance vector: ", result)




# Test cases
# v1 = Vector(3)
# v2 = Vector(3)
# v1.set(0, 1)
# v1.set(1, 2)
# v1.set(2, 10)
# v2.set(0, 1)
# v2.set(1, 5)


# print(v1 + v2)  # Vector(2, 7, 10)

# v3 = v1.shift(1)
# print(v3)  # Vector(0, 1, 2)