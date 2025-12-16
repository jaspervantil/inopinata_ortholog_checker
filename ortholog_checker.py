import pandas as pd
import re
from Bio import SeqIO
from Bio import Blast
import subprocess

#open database as df
database=pd.read_excel(r"C:\Ortholog_Checker\elegans_inopinata orthologs.xlsx")
elegans_proteome=r"C:\Ortholog_Checker\Blast_Temp\elegans_proteome.fasta"
inopinata_proteome=r"C:\Ortholog_Checker\Blast_Temp\inopinata_proteome.fa"


#ask input and make query fasta
POI = input("enter gene")
if POI in database["elegans gene"].values:
    with open(elegans_proteome) as e_handle:
        for e_record in SeqIO.parse(e_handle, "fasta"):
            if POI in e_record.description:
                match = re.search(r'GN=([^\s]+)',e_record.description)
                GN=match.group(1) if match else None
                e_seq=e_record.seq
    with open(inopinata_proteome) as i_handle:
        for i_record in SeqIO.parse(i_handle, "fasta"):
            ortholog_id=database.loc[database['elegans gene'] == POI, 'inopinata ortholog id'].iloc[0]
            ortholog_id = ortholog_id.strip('"').strip("'")
            if i_record.id.startswith(ortholog_id):
                i_seq=i_record.seq
    print("Existing Entry:", GN)
    print("Score:", database.loc[database['elegans gene']==POI,'Score'], "%Identity:", database.loc[database['elegans gene']==POI,'%Identity'])
    print('elegans', POI, 'sequence:',e_seq)
    print('inopinata ortholog sequence:',i_seq)

else:
    with open(elegans_proteome) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if POI in record.description:
                e_id=record.id
                e_seq=record.seq
                SeqIO.write(record, r'C:\Ortholog_Checker\Blast_Temp\query_tmp.fasta', 'fasta')
    cmd=r'blastp -query C:\Ortholog_Checker\Blast_Temp\query_tmp.fasta -db C:\Ortholog_Checker\Blast_Temp\inopinata_proteome_db -out C:\Ortholog_Checker\Blast_Temp\output.xml -outfmt 5'
    subprocess.run(cmd,shell=True)
    blast_record = Blast.read(r"C:\Ortholog_Checker\Blast_Temp\output.xml")
    hit=blast_record[0]
    alignment=hit[0]
    o_id=hit.target.description.split(' ', 1)[0]
    with open(inopinata_proteome) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id.startswith(o_id):
                o_seq=record.seq
    o_score=alignment.annotations["bit score"]
    o_perid = 100 * alignment.annotations["identity"] / len(hit.target.seq)
    out_df=pd.DataFrame([{'elegans gene':POI,
                         'elegans gene id':e_id,
                         'inopinata ortholog id':o_id,
                         'Score':int(o_score),
                         '%Identity':int(o_perid)}])
    new_db=pd.concat([database,out_df], ignore_index=True)
    new_db.to_excel(r"C:\Ortholog_Checker\elegans_inopinata orthologs.xlsx", index=False)
    print("New Entry:", POI)
    print("Score:", int(o_score), "%Identity:", int(o_perid))
    print("elegans", POI, "sequence:", e_seq)
    print("inopinata ortholog sequence", o_seq)