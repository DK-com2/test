import os
import shutil
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

folder_path = '/home/dkcom/デスクトップ/primer_design/Sympetrum_Nucleotide_number/'

for file_name in os.listdir(folder_path):
    if file_name.endswith('.seq'):
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, 'r') as f:
            ids = []
            for line in f:
                ids.append(line.strip())
        name = file_name.replace('_COI.seq','')
        print(name)
        print(ids)

        # EntrezのメールアドレスとAPIキーを設定する

        Entrez.email = "y200561@mail.ryukoku.ac.jp"
        Entrez.api_key =  "c6deea20447a21c35cc36f7b8699bf343708"

        # Entrez.efetchを使って配列データを取得する
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))

        # 配列データを保存する
        SeqIO.write(records, "test.fasta", "fasta")

        # アライメントする配列を指定する
        clustalw_cline = ClustalwCommandline("clustalw", infile="test.fasta")

        # アライメントを実行する
        stdout, stderr = clustalw_cline()

        # アライメント結果を読み込む
        alignment = AlignIO.read("test.aln", "clustal")


        # .fastaファイルに書き出す
        new_records = SeqIO.write(alignment, name+".fasta", "fasta")

        #不要なファイルの削除
        os.remove("test.aln")
        os.remove("test.dnd")
        os.remove("test.fasta")


# フォルダのパス
folder_path = '/home/dkcom/デスクトップ/primer_design/alignment_python/'
# "_"で区切って、最初の部分を取り出す
first_part = name.split("_")[0]

print(first_part)  # Output: "Sympetrum"

# フォルダが存在しない場合に作成
if not os.path.exists(first_part):
    os.makedirs(first_part)


    # Sympetrum_~.fastaをSympetrumフォルダに移動
    # 現在のフォルダ内にあるSympetrum_~~.fastaファイルをすべて取得
    for filename in os.listdir():
        if filename.startswith(first_part) and filename.endswith(".fasta"):
            # Sympetrumフォルダに移動
            shutil.move(filename, first_part)

else:
    for filename in os.listdir(folder_path):
        if filename.startswith(first_part+"_") and filename.endswith(".fasta"):
            os.remove(os.path.join(filename))
