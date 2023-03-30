# from Bio.Phylo.PAML import codeml
# import os

# # Parse the codeml output file
# results = codeml.read("SCYL3_codeml_test.out")
# with open('SCYL3.test', 'w') as fp:

#     d=results.readlines()
#     results.seek(0)
#     for i in d:
#         if i == "lnL(...":
#             fp.write(i)




# # Open the input and output files
# with open("SCYL3_codeml_test.out", 'r') as input_file, open('lnL_output.txt', 'w') as output_file:
    
#     # Loop through the lines in the input file
#     for line in input_file:
        
#         # If the line starts with "lnL", extract the value and write it to the output file
#         if line.startswith('lnL'):
#             lnL_value = float(line.split()[4])
#             output_file.write(str(lnL_value) + '\n')




with open('SCYL3_codeml_test.out') as f:
    model_num = ''
    for line in f:
        if line.startswith('Model'):
            model_num = line.split()[1].strip(':')
        elif line.startswith('lnL'):
            lnL_value = float(line.split()[4])
            with open('output2.txt', 'a') as out_file:
                out_file.write(f'Model {model_num}: lnL = {lnL_value}\n')
