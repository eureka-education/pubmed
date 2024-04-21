# 必要なライブラリのインポート
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import streamlit as st

# Entrezにメールアドレスを登録します（必須）
Entrez.email = "your_email@example.com"

# PubMedから論文の要約を取得する関数
def fetch_summary_from_pubmed(doi):
    try:
        # DOIを使用してPubMedから論文のHTMLを取得
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={doi}"
        response = requests.get(pubmed_url)
        if response.status_code == 200:
            # BeautifulSoupを使用してHTMLを解析
            soup = BeautifulSoup(response.content, "html.parser")
            # 要約を含む<div>要素を特定
            abstract_div = soup.find("div", class_="abstract-content")
            # 要約のテキストを抽出
            if abstract_div:
                abstract_text = abstract_div.text.strip()
                return abstract_text
        return "要約がありません"
    except Exception as e:
        return f"エラーが発生しました: {e}"

# PubMedから論文を検索し、論文の情報と要約を取得する関数
def search_and_fetch_pubmed_articles(query, max_results=5):
    try:
        # 検索クエリを使用してPubMedから論文を検索（最新順に最大5件）
        handle = Entrez.esearch(db="pubmed", term=query, sort="most recent", retmax=max_results)
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        articles = []
        for pubmed_id in id_list:
            summary_handle = Entrez.esummary(db="pubmed", id=pubmed_id)
            summary_record = Entrez.read(summary_handle)[0]
            
            # タイトルを取得
            title = summary_record.get("Title", "タイトルがありません")
            
            # DOIを取得
            doi = summary_record.get("DOI", "DOIがありません")
            
            # 著者を取得
            authors = ', '.join(summary_record.get("AuthorList", []))
            
            # DOIを使用して要約を取得
            abstract = fetch_summary_from_pubmed(doi)
            
            article = {
                "title": title,
                "authors": authors,
                "abstract": abstract,
                "doi": doi
            }
            articles.append(article)
        
        return articles
    except Exception as e:
        return []

# Streamlitアプリのタイトルを設定
st.title("PubMed論文検索")

# 検索クエリを入力するテキストボックスを表示
query = st.text_input("検索クエリを入力してください")

# 検索ボタンがクリックされたら、PubMed APIを使用して論文を検索し、要約を取得して表示
if st.button("検索"):
    if query:
        articles = search_and_fetch_pubmed_articles(query, max_results=5)
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
    else:
        st.warning("検索クエリを入力してください。")


