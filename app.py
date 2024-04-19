from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import streamlit as st

# Entrezにメールアドレスを登録します（必須）
Entrez.email = "l2hhok8g40@sute.jp"

# PubMedから論文を検索する関数
def search_pubmed(query):
    try:
        # 検索クエリを使用してPubMedから論文を検索
        handle = Entrez.esearch(db="pubmed", term=query)
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        articles = []
        for pubmed_id in id_list:
            summary_handle = Entrez.esummary(db="pubmed", id=pubmed_id)
            summary_record = Entrez.read(summary_handle)[0]
            
            # タイトルと要約を取得
            title = summary_record.get("Title", "タイトルがありません")
            # 著者を取得
            authors = ', '.join(summary_record.get("AuthorList", []))
            # DOIを取得
            doi = summary_record.get("DOI", "DOIがありません")

            # 要約を取得するための関数を呼び出し
            abstract = fetch_summary_from_doi(doi)
            
            article = {
                "title": title,
                "authors": authors,
                "abstract": abstract,
                "doi": doi
            }
            articles.append(article)
        
        return articles
    except Exception as e:
        st.error(f"エラーが発生しました: {e}")

# 論文のDOIから要約を取得する関数
def fetch_summary_from_doi(doi):
    # DOI解決サービスを使用して論文の情報を取得
    crossref_api_url = f"https://api.crossref.org/works/{doi}"
    response = requests.get(crossref_api_url)
    if response.status_code == 200:
        data = response.json()
        if 'abstract' in data['message']:
            abstract = data['message']['abstract']
            return abstract
    return "要約がありません"

# Streamlitアプリのタイトルを設定
st.title("Pubmed論文検索")

# 検索クエリを入力するテキストボックスを表示
query = st.text_input("検索クエリを入力してください")

# 検索ボタンがクリックされたら、PubmedAPIを使用して論文を検索
if st.button("検索"):
    # PubmedAPIを使用して論文を検索し、結果を表示
    articles = search_pubmed(query)
    
    if articles:
        st.write(f"検索結果: {len(articles)} 件")
        
        for article in articles:
            st.write("-------")
            st.write(f"タイトル: {article['title']}")
            st.write(f"著者: {article['authors']}")
            st.write(f"要約: {article['abstract']}")
            st.write(f"DOI: {article['doi']}")
    else:
        st.warning("検索結果が見つかりませんでした。別のクエリをお試しください。")
