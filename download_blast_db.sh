#!/bin/bash

# スクリプトの途中でエラーが発生したら処理を停止します
set -e

# --- 関数定義 ---
# スクリプトの使い方を表示する関数です
usage() {
    echo "NCBIからBLASTデータベースをダウンロードし、その場で展開します。"
    echo ""
    echo "Usage: $0 <database_name> [output_directory]"
    echo ""
    echo "Arguments:"
    echo "  <database_name>    : ダウンロードしたいデータベース名 (例: swissprot, nr)"
    echo "  [output_directory] : (任意) 保存先のディレクトリ。省略した場合はカレントディレクトリに保存します。"
    echo ""
    echo "Example (カレントディレクトリに保存):"
    echo "  $0 swissprot"
    echo ""
    echo "Example (指定ディレクトリに保存):"
    echo "  $0 nr /home/masakazu/blast/db"
    exit 1
}

# --- メイン処理 ---

# コマンドライン引数（データベース名）が指定されているかチェックします
if [ -z "$1" ]; then
    echo "エラー: データベース名が指定されていません。"
    usage
fi

# --- 変数の設定 ---
DB_NAME="$1"
# 第2引数が指定されていればそれを保存先に、なければカレントディレクトリ「.」を保存先にします
OUTPUT_DIR="${2:-.}"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/blast/db/"

echo "----------------------------------------------------"
echo " BLASTデータベースのセットアップを開始します"
echo " データベース名: ${DB_NAME}"
echo " 保存先ディレクトリ: $(realpath "${OUTPUT_DIR}")"
echo "----------------------------------------------------"

# 保存先ディレクトリを作成し、そこに移動します
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

# --- ファイルリストの取得 ---
echo "ステップ1: ダウンロード対象のファイルリストを取得します..."
# .tar.gz ファイルのリストだけを取得します
TAR_FILE_LIST=$(curl -sL "${BASE_URL}" | grep -oP 'href="('"${DB_NAME}"'[^"]*\.tar\.gz)"' | sed -e 's/href="//' -e 's/"//' || true)

# ダウンロード対象が見つからない場合はエラーで終了します
if [ -z "$TAR_FILE_LIST" ]; then
    echo ""
    echo "エラー: ${DB_NAME} という名前のデータベースファイルが見つかりませんでした。"
    echo "データベース名が正しいか、再度確認してみてください。"
    exit 1
fi

# --- ダウンロード、検証、解凍のループ処理 ---
echo "ステップ2: ダウンロード、検証、解凍をファイルごとに実行します..."
for tar_file in $TAR_FILE_LIST; do
    md5_file="${tar_file}.md5"
    echo ""
    echo ">>> 処理中のファイル: ${tar_file}"

    # 1. 圧縮ファイル (.tar.gz) をダウンロード
    echo "  (1/4) 圧縮ファイルをダウンロード中..."
    wget -q --show-progress -c "${BASE_URL}${tar_file}"

    # 2. チェックサムファイル (.md5) をダウンロード (存在すれば)
    echo "  (2/4) チェックサムファイルをダウンロード中..."
    # 404エラーでスクリプトが停止しないように `|| true` をつけます
    wget -q -c "${BASE_URL}${md5_file}" || true

    # 3. チェックサムを検証 (md5ファイルが存在すれば)
    echo -n "  (3/4) ファイルの整合性をチェック中... "
    if [ -f "$md5_file" ]; then
        if md5sum -c --quiet "$md5_file"; then
            echo "OK"
        else
            echo "警告: 整合性チェックに失敗しました！ (${tar_file}) スクリプトを停止します。"
            exit 1 # 壊れたファイルの可能性があるので、ここで停止する方が安全です
        fi
    else
        echo "スキップ (チェックサムファイルなし)"
    fi

    # 4. ファイルを解凍
    echo "  (4/4) ファイルを解凍中..."
    tar -xvzf "$tar_file"

    # 5. 後片付け
    echo "        -> 処理済みの圧縮ファイルとチェックサムファイルを削除します。"
    rm -f "$tar_file"
    rm -f "$md5_file"
done

echo ""
echo "🎉🎉🎉 全ての処理が完了しました！ 🎉🎉🎉"
echo "データベースファイルは $(pwd) に展開されています。"
