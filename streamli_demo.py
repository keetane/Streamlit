# https://camp.trainocate.co.jp/magazine/streamlit-web/
import streamlit as st

st.title("title")
st.write('write ' + "write")
st.markdown("# Head1")
st.markdown("## Head2")

# st.checkbox("チェックボックス") #引数に入れることでboolを返す
# st.button("ボタン") #引数に入れるとboolで返す
# st.selectbox("メニューリスト", ("選択肢1", "選択肢2", "選択肢3")) #第一引数：リスト名、第二引数：選択肢
# st.multiselect("メニューリスト（複数選択可）", ("選択肢1", "選択肢2", "選択肢3")) #第一引数：リスト名、第二引数：選択肢、複数選択可
# st.radio("ラジオボタン", ("選択肢1", "選択肢2", "選択肢3")) #第一引数：リスト名（選択肢群の上に表示）、第二引数：選択肢
# st.text_input("文字入力欄") #引数に入力内容を渡せる
# st.text_area("テキストエリア") #引数に入力内容を渡せる


# check = st.checkbox("チェックボックス") #引数に入れることでboolを返す

# if check:
#    st.button("ボタン") #引数に入れるとboolで返す
#    st.selectbox("メニューリスト", ("選択肢1", "選択肢2", "選択肢3")) #第一引数：リスト名、第二引数：選択肢
#    st.multiselect("メニューリスト（複数選択可）", ("選択肢1", "選択肢2", "選択肢3")) #第一引数：リスト名、第二引数：選択肢、複数選択可
#    st.radio("ラジオボタン", ("選択肢1", "選択肢2", "選択肢3")) #第一引数：リスト名（選択肢群の上に表示）、第二引数：選択肢
#    st.text_input("文字入力欄") #引数に入力内容を渡せる
#    st.text_area("テキストエリア")

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
# 環状化合物をきれいに描画 https://future-chem.com/rdkit-coordgen/
rdDepictor.SetPreferCoordGen(True)

# %%
import streamlit as st
check = st.checkbox("チェックボックス") #引数に入れることでboolを返す
button = st.button("ボタン", 'test') #引数に入れるとboolで返す
if button:
    check
# %%
st.divider()
select = ["test", "選択肢2", "選択肢3"]
st.selectbox("メニューリスト", (select)) #第一引数：リスト名、第二引数：選択肢
st.multiselect("メニューリスト（複数選択可）", (select)) #第一引数：リスト名、第二引数：選択肢、複数選択可
st.radio("ラジオボタン", (select)) #第一引数：リスト名（選択肢群の上に表示）、第二引数：選択肢
# 以下をサイドバーに表示
st.sidebar.text_input("文字入力欄") #引数に入力内容を渡せる
text = st.text_area("input SMILES ")

st.write(text,)
if text:
    # st.write(text)
    st.image(Draw.MolToImage(Chem.MolFromSmiles(text)), width=400)


