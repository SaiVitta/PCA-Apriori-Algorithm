"""
Description: Get number of frequent sets for all the lengths of sets at every step of Apriori Algorithm
Group Members: Mitali Bhiwande | Sumedh Ambokar | Tejasvi Sankhe
"""

import numpy as np
import sys


# Function to read the data file. Every gene is represented using proper naming convention. {GNumber_gene}
def read_data(file_name):
    with open(file_name) as input_data:
        genes=[line.replace("\n","").split('\t') for line in input_data]
    dataset_list=[]
    for sample in np.array(genes):
        i = 1
        data_set = set()
        for gene in sample:
            data_set.add("G" + str(i) + "_" + gene.replace(" ", "_"))
            i = i + 1
        dataset_list.append(data_set)
    return dataset_list


# Function to get support count(Number of occurences) of a set from the given data set
def get_support_count(input_list, dataset_list):
    count = 0
    for data_set in dataset_list:
        if set(input_list).issubset(data_set):
            count = count + 1
    return count


# Function used to generate new sets of size k + 1 from the previous frequent sets which are of size k.
def generate_sets(size_of_set, previous_freq_dataset_list, dataset_list, support_count):
    i = 0
    list_of_sets = []
    flag = False
    if size_of_set == 1:
        dic={}
        for sample in dataset_list:
            for gene in sample:
                if gene in dic:
                    dic.update({gene:dic.get(gene) +1})
                else: 
                    dic.update({gene:1})
        for key in list(dic.keys()):
            if dic.get(key) < support_count:
                del dic[key]      
        return list(dic)
    else:
        for first_dataset in previous_freq_dataset_list[i:]:
            for second_dataset in previous_freq_dataset_list[i+1:]:
                if size_of_set == 2:
                    pair = []
                    pair.append(first_dataset)
                    pair.append(second_dataset)
                    if get_support_count(pair, dataset_list) >= support_count:
                        list_of_sets.append(pair)
                else:
                    for ele in range(size_of_set - 2):
                        if first_dataset[ele] != second_dataset[ele]:
                            flag = True
                            break
                    if flag:
                        flag = False
                        break
                    else:
                        dataset = []
                        dataset = dataset+first_dataset
                        dataset.append(second_dataset[len(second_dataset) - 1])
                        sc = get_support_count(dataset, dataset_list)
                        if sc >= support_count:
                            list_of_sets.append(dataset)
            i = i + 1
    return list_of_sets 


# Function to calculate the number of frequent set of every length given a support count percentage
def calulate_frequent_set_count(support_count_percent, dataset_list):
    frequency_dic = {}
    support_count = int(support_count_percent*len(dataset_list)/100)
    freq_set_length_count = 1
    previous_freq_dataset_list = []
    set_size = 1
    while freq_set_length_count > 0:
        previous_freq_dataset_list = generate_sets(set_size, previous_freq_dataset_list, dataset_list, support_count)
        freq_set_length_count = len(previous_freq_dataset_list)
        if freq_set_length_count > 0:
            frequency_dic.update({set_size : freq_set_length_count})
        set_size += 1
    return frequency_dic


# Reads parameters from the command line inputs.
if len(sys.argv) > 1:	
	dataset_list = read_data(sys.argv[1])
	support = int(sys.argv[2])
	length_dic = calulate_frequent_set_count(support, dataset_list)
	for length in length_dic:
		print("number of length",length ,"frequent itemsets:", length_dic[length])
else:
	print("Not sufficient Inputs")

