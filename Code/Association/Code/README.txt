###################### Instructions to execute code ######################


### For Part A:
- The code file takes in 2 parameters 1. Data file path 2. Support Count Percentage
- Below is a sample execution command:
>> python CSE601_Project_1_Part_2A.py associationruletestdata.txt 50


### For Part B:
- The code file takes varied number of inputs depending on the template used.
- Below is list of common parameters for all the 3 templates:
	1. Data file path
	2. Support Count Percentage
	3. Confidence Percentage
	4. Template type. i.e. template1 or template2 or template3
- All other parameters are passed after these four in the sequence provided by the template.
- Genes inputs need to be passed as string seperated by commas. For example: ["G1_Up","G1_Down"] -> G1_Up,G1_Down
- Format of Disease has been changed to maintain consistency in data. Disease is treated as G101_disease_name. For example: "ALL" --> "G101_ALL"
- If there exists a space in the disease name, for example: "Breast Cancer", then replace space by "_" while passing as parameter. Thus input for "Breast_Cancer" --> "G101_Breast_Cancer"
- Below are sample execution commands for each template:
1. Template 1: template1("RULE","ANY",["G1_Up","G1_Down"])
>> python CSE601_Project_1_Part_2A.py associationruletestdata.txt 50 70 template1 RULE ANY G1_Up,G1_Down
>> python CSE601_Project_1_Part_2A.py associationruletestdata.txt 50 70 template1 RULE ANY G1_Up,G101_Breast_Cancer
2. Template 2: template2("RULE",3)
>> python CSE601_Project_1_Part_2A.py associationruletestdata.txt 50 70 template2 RULE 3
3. Template 3: template3("1or2","RULE","ANY",["G1_Up"],"RULE",3)
>> python CSE601_Project_1_Part_2A.py associationruletestdata.txt 50 70 template3 1or2 RULE ANY G1_Up RULE 3