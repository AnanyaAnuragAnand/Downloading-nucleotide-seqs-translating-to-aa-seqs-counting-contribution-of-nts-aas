import pandas
import openpyxl
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

Entrez.email = 'ananyaanurag12@gmail.com'

raw_data = pandas.read_excel("./acc_new.xlsx", "Sheet1")

new_list = []
for name in raw_data["id"]:
    input_id = name.strip()
    new_list.append(input_id)
    handle = Entrez.efetch(db="nucleotide", id=new_list, rettype="fasta")
    record = SeqIO.parse(handle, "fasta")
    outputname = "./seq_string.fasta"
    SeqIO.write(record, outputname, "fasta")

seq_data = open(file="./seq_string.fasta", mode="r")
# for line in seq_data.read().split("\n"):
#     seqs += line
seqs = seq_data.read().split("\n")
s_list = []
append_data = ''
for rec in seqs:
    if rec == '':
        pass
    if '>' in rec:
        if append_data:
            s_list.append(append_data)
            append_data = ''
        s_list.append(rec)

    else:
        append_data += rec
s_list.append(append_data)

keys, vals = [], []
for index in range(0, len(s_list), 2):
    keys.append(s_list[index])

for index in range(1, len(s_list), 2):
    vals.append(s_list[index])

# print("!!", len(keys), len(vals), len(s_list))
# print("!!", s_list)

final, data_frame_list, data_frame_list_tb = [], [], []
for x in range(0, 154):
    read_list, read_list_tb = [], []
    body = Seq(vals[x])
    A = body.count("A")
    T = body.count("T")
    G = body.count("G")
    C = body.count("C")
    nt_count = open("./nt_count.txt", "a")

    xl = open("C:/Users/anany/OneDrive/Desktop/Computational_Biology/acc_ids.txt", "r")
    xl_file = xl.read().split("\n")

    nt_count.write(f"{A}, {T}, {G}, {C} \n")
    read_list.append(xl_file[x])
    read_list.append(A)
    read_list.append(T)
    read_list.append(G)
    read_list.append(C)
    data_frame_list.append(read_list)
    nt_count.close()


    header = keys[x]
    b = body.translate()


    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "O", "S", "U", "T", "W", "Y", "V", "B", "Z", "X", "J"]
    aa_count = open("./aa_count.txt", "a")
    read_list_tb.append(xl_file[x])
    for a in aa:
        a_count = b.count(a)
        aa_count.write(f"{a_count}, ")
        read_list_tb.append(a_count)
    data_frame_list_tb.append(read_list_tb)
    aa_count.write(f"\n ")
    c = f"{header} \n {b}"
    final.append(c)

# writer = pandas.ExcelWriter('ntdf.xlsx', engine='xlsxwriter')
# writer.save()
df = pandas.DataFrame(data_frame_list, columns=['id', 'A', 'T', 'G', 'C'])
df.to_excel('ntdf.xlsx', index=False)

df = pandas.DataFrame(data_frame_list_tb, columns=["Id", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "O", "S", "U", "T", "W", "Y", "V", "B", "Z", "X", "J"])
df.to_excel('aadf.xlsx', index=False)


new_seq = pandas.DataFrame(final)
a = new_seq.squeeze()
# print(a[0])
f = open(f"./seqt.fasta", "w")
f.write(f"{a}")
f.close()

######################################################################################

