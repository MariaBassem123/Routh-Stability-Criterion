import numpy as np
from sympy import *

# -------------------- Examples --------------------
# There are more test cases in the "test cases.txt" file.
# s^5+s^4+10s^3+72s^2+152s+240 ==> in the project (unstable) [2 roots in RHS]
# 2s^5+3s^4+2s^3+3s^2+2s+1 ==> epsilon --> unstable [2 root in RHS]
# s^5+2s^4+24s^3+48s^2-25s-50 ==> whole row = 0 (unstable) [1 in RHS]
# s^5+7s^4+6s^3+42s^2+8s+56 ==> whole row = 0 (marginally stable) [1 LHS, 4 jw]
# s^4+2s^3+6s^2+4s+1 --> stable system

# -------------------- Global Variables --------------------
coefficients = []
powers_of_s = []
exponents = []
routh_array = []
epsilon = 1e-10
zero_row = 0 # zero means no row equal zero
imj_count = 0

# -------------------- function definitions --------------------
def print_routh_table(routh_array,powers_of_s):
    for i in range(0,len(routh_array)):
        print(powers_of_s[i],"| ",end="")
        for j in range(0,len(routh_array[i])):
            print(routh_array[i][j]," ",end="")
        print("\n")

def parse_input(input):
    if input.find('+') != -1:
        parsed = input.replace(' ','').replace('+', ' ').replace('-', ' -').replace('**', '^').replace('*', '').split()
        # print("parsed input = ",parsed)
        coeffs = []
        powers_of_s = []
        exponents = []
        for str in parsed:
            coeff = ''
            s_power = ''
            exponentFlag = False # not an s^ sth
            for i in str:
                if i != 's' and not exponentFlag:
                    coeff += i
                else:
                    s_power += i
                    exponentFlag = True
            if coeff == '':
                coeff = '1'
            if s_power == '':
                s_power = 's^0'
            elif s_power == 's':
                s_power = 's^1'
            coeffs.append(coeff)
            powers_of_s.append(s_power)
        for i in powers_of_s:
            i = i.replace('s^','')
            exponents.append(i)
    return coeffs,powers_of_s,exponents

def initialize_routh_array(coeffs):
    first_row = []
    second_row = []
    i = 0
    for col in coeffs:
        if i%2 == 0:
            first_row.append(col)
        else:
            second_row.append(col)
        i += 1

    while len(second_row) < len(first_row):
        second_row.append(0)
    first_row = [float(i) for i in first_row]
    second_row = [float(i) for i in second_row]
    return first_row, second_row

def fix_zero_row(row,previous):
    aux_order = 0
    # The whole row equals zero. So we replace it with the derivative of the previous row
    global zero_row # refers to the global variable
    zero_row += 1
    index = len(routh_array)-1
    cur_pow = exponents[index]
    equation = ''
    for x in range(0,len(previous)):
        equation += str(previous[x]) + '*s**' + str(cur_pow)
        cur_pow = int(cur_pow) - 2
        if not (x == len(previous)-1):
            if x < len(previous) and np.sign(previous[x+1]) == -1:
                continue
            equation += '+'
        if cur_pow < 0:
            continue
    s = symbols('s')
    df = diff(equation,s)
    derivative = parse_input(str(df))
    exp = [int(i) for i in derivative[2]]
    aux_order = max(exp) + 1
    row = derivative[0]
    while len(row) < len(previous):
        row.append(0)
    row = [float(i) for i in row]
    return row,aux_order

def eval_next_row(previous, before_previous):
    # given the last 2 rows in the table, we can compute the next row
    # previous row doesn't have all zeros (Normal Computation)
    # before_previous ==> c   d
    # previous ==> a    b
    row = []
    a = previous[0]
    c = before_previous[0]
    if a == 0.0:
        #Replace it with very small number epsilon because we can't divide by 0
        a = epsilon
    for col in range(0, len(previous)):
        val = 0.0
        if (col + 1) < len(previous):
            d = before_previous[col+1]
            b = previous[col+1]
            val = (a*d-b*c)/a
        row.append(val)
        count = 0
    for i in row:
        if i == 0: 
            count += 1
    global imj_count
    # check if a row has all elements = zero 
    if count == len(row):
        row,imj_count = fix_zero_row(row,previous)
    return row,imj_count

def isStable(routh_array):
    first_column = []
    num_rows = len(routh_array)
    for x in range(0,num_rows):
        first_column.append(routh_array[x][0])
    
    flag = 'yes' # stability flag
    c = 0
    count = 0
    sign = np.sign(first_column[0])
    for i in range(1,len(first_column)):
        global zero_row
        if (np.sign(first_column[i]) == sign and zero_row > 0) or (np.sign(first_column[i]) == 0):
            flag = 'could be marginal'
        elif not(np.sign(first_column[i]) == sign):
            flag = 'no'
            c += 1
        elif np.sign(first_column[i]) == sign and c > 0:
            count =+ c
            c = 0
        if not(np.sign(first_column[i]) == 0):
            sign = np.sign(first_column[i])
    if not(c == 0):
        count =+ c
        c = 0
    return first_column,flag,count


# -------------------- Beginning of the program --------------------   
input = input('Enter characteristic equation:\n')
coefficients,powers_of_s,exponents = parse_input(input)
first_row, second_row = initialize_routh_array(coeffs= coefficients)
routh_array.append(first_row)
routh_array.append(second_row)
system_order = max(exponents)

if system_order == 1:
    # no need to compute other rows, only the first and second rows
    # checking the stability of the system
    first_column = isStable(routh_array, coefficients)

else:
    for i in range(0,len(powers_of_s)-2):    
        new_row,imj_count = eval_next_row(routh_array[-1], routh_array[-2])
        # print("new_row = ",new_row)
        routh_array.append(new_row)
first_column,flag,RHS_count = isStable(routh_array)

print("=========================================")
print("Answer: ")
print("Routh Array = ")
print_routh_table(routh_array,powers_of_s)
print("System State: ", end="")
if flag == 'yes':
    print("Stable")
elif flag == 'no':
    print("Unstable")
else: # Marginal
    print("Marginally stable")

print("Number of roots in RHS: ", RHS_count)
print("Number of roots in LHS: ", int(system_order) - (RHS_count+imj_count))
print("Number of roots on jw axis: ",imj_count)