import os
import shutil
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from IPython.display import clear_output
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from selenium import webdriver
from pprint import pprint





#アクセッション番号リストのフォルダパス
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


# alignment_pythonフォルダのパス
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

#########################################################################################################################



#相対パスを指定
folder_alignment_path = first_part


#フォルダ内のファイルを一つずつ読み込む
for file_name in os.listdir(folder_alignment_path):
    if file_name.endswith('.fasta'):
        file_path = os.path.join(folder_alignment_path, file_name)
        alignment = AlignIO.read(file_path, 'fasta')
        #一連の処理の進捗を確認するために表示
        name = file_name.replace('.fasta','')
        print(name)
        time.sleep(0.5)
        


#代表配列を保存するための空のリストを作成する
        consensus = []
        for i in range(alignment.get_alignment_length()):
        # 配置の列を取得する
            column = alignment[:, i]
        # 列の A、C、G、T、- の数を数える
            count = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
            for j in column:
                if j in count:
                    count[j] += 1
        # カウントが最も高い値を返す
            max_base = max(count, key=count.get)
        # 配列に塩基を追加する
            consensus.append(max_base)


        #代表配列の SeqRecord オブジェクトを作成する
        consensus_seq = Seq("".join(consensus))
        consensus_record = SeqRecord(consensus_seq)

        #一連の処理の進捗を確認するために表示
        print(consensus_record.format("fasta"))
        time.sleep(0.5)
        clear_output(wait=True)

        #代表配列を保存する
        consensus_record.id = name
        consensus_record.description = 'consensus sequence'
        SeqIO.write(consensus_record, name + '_consensus.fasta', 'fasta')


# "_"で区切って、最初の部分を取り出す
first_part = name.split("_")[0]

# フォルダが存在しない場合に作成
if not os.path.exists(first_part+"_consensus"):
    os.makedirs(first_part+"_consensus")
    # フォルダに移動
    shutil.move(name + '_consensus.fasta', first_part+"_consensus")

# フォルダが存在する場合に移動

    for filename in os.listdir():
        if filename.startswith(first_part+"_") and filename.endswith("_consensus.fasta"):
        # フォルダに移動
            shutil.move(filename, first_part+"_consensus")
else:
    for filename in os.listdir(folder_path):
        if filename.startswith(first_part+"_") and filename.endswith("_consensus.fasta"):
            os.remove(os.path.join(filename))



#########################################################################################################################



#taget種のファイルパス
def select_file(path):
    files = os.listdir(path)
    for i, file in enumerate(files, start=1):
        print(f"{i}. {file}")
    choice = None
    while choice is None:
        try:
            choice = int(input("Enter the number of the file you want to select: "))
            if choice not in range(1, len(files) + 1):
                raise ValueError
        except ValueError:
            print("Invalid option, please enter a number between 1 and", len(files))
            choice = None
    selected_file = files[choice - 1]
    return os.path.join(path, selected_file)

target_path = first_part+"_consensus"

#黄色文字で表示
print("\033[33m{}\033[0m".format(f"ターゲット種の配列データを選択してください"))
selected_file_path = select_file(target_path)
print("\033[33m{}\033[0m".format(f"Selected file: {selected_file_path}"))

target = AlignIO.read(selected_file_path, 'fasta')
print(target.format('fasta'))


#クロームの立ち上げ
driver=webdriver.Chrome()

#URLの指定
driver.get("https://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plus.cgi")
time.sleep(1)

#検索ボックスの指定
search_box=driver.find_element_by_xpath("/html/body/div/div[2]/form/div[5]/div[3]/table[1]/tbody/tr[2]/td/textarea")

#検索ボックスに文字を入力
search_box.send_keys(target.format('fasta'))

#検索ボタンの指定
search_button=driver.find_element_by_xpath("/html/body/div/div[2]/form/div[2]/table/tbody/tr/td[3]/table/tbody/tr/td[1]/input")

#検索ボタンをクリック
search_button.click()
time.sleep(0.5)

#ブラウザのウインドウ高を取得する
win_height = driver.execute_script("return window.innerHeight")

#スクロール開始位置の初期値（ページの先頭からスクロールを開始する）
last_top = 1

#ページの最下部までスクロールする無限ループ
while True:

  #スクロール前のページの高さを取得
  last_height = driver.execute_script("return document.body.scrollHeight")
  
  #スクロール開始位置を設定
  top = last_top

  #ページ最下部まで、徐々にスクロールしていく
  while top < last_height:
    top += int(win_height * 0.1)
    driver.execute_script("window.scrollTo(0, %d)" % top)
    time.sleep(0.05)

  #１秒待って、スクロール後のページの高さを取得する
  time.sleep(1)
  new_last_height = driver.execute_script("return document.body.scrollHeight")

  #スクロール前後でページの高さに変化がなくなったら無限スクロール終了とみなしてループを抜ける
  if last_height == new_last_height:
    break

  #次のループのスクロール開始位置を設定
  last_top = last_height

# primer dict
primers = {}
i = 0
x = 0
while True:
    try:
        forward_primer = driver.find_element_by_id(f"PRIMER_{i}_SEQUENCE").get_attribute("value")
        i += 1
        x += 1
        reverse_primer = driver.find_element_by_id(f"PRIMER_{i}_SEQUENCE").get_attribute("value")
        primers[f"primer{x}"] = {"forward primer": forward_primer, "reverse primer": reverse_primer}
        i += 1
    except:
        break

clear_output(wait=True)
pprint(primers)

#クロームの終了
driver.quit()

#########################################################################################################################

#クロームの立ち上げ
driver=webdriver.Chrome()

#URLの指定
driver.get("https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?")
time.sleep(3)

for primer_name, primer_values in primers.items():
    
    #forward primerの入力
    search_box=driver.find_element_by_xpath("/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[2]/div[1]/input")
    search_box.send_keys(primer_values['forward primer'])
    time.sleep(0.2)
    #reverse primerの入力
    search_box=driver.find_element_by_xpath("/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[2]/div[2]/input")
    search_box.send_keys(primer_values['reverse primer'])

    #少しスクロール
    time.sleep(0.3)
    driver.execute_script("window.scrollTo(0, 600)")

    #data baseの指定
    database_select = driver.find_element_by_xpath("/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div[3]/span/select")
    database_select.send_keys("nr")
    time.sleep(0.3)
    #Organismを空にする
    search_button=driver.find_element_by_xpath("/html/body/div[1]/div[2]/div[2]/form/div[2]/fieldset[4]/div[5]/div/div/input[1]")
    search_button.clear()

    time.sleep(0.3)
    #Get primerをクリック
    search_button=driver.find_element_by_xpath("/html/body/div[1]/div[2]/div[2]/form/div[3]/div[1]/input")
    search_button.click()

    # 新しいタブを開く
    driver.execute_script("window.open('https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?');")
    # 新しいタブに切り替え
    driver.switch_to.window(driver.window_handles[-1])
    time.sleep(30)



