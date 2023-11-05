import streamlit as st

# ページのタイトル設定
st.set_page_config(
    page_title="Wedding Quiz",
)

# セッション情報の初期化
if "page_id" not in st.session_state:
    st.session_state.page_id = "main"
    st.session_state.answers = []

# 各種メニューの非表示設定
hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden; }
        footer {visibility: hidden;}
        </style>
        """
st.markdown(hide_menu_style, unsafe_allow_html=True)

# 最初のページ
def main():
    st.markdown(
        "<h1 style='text-align: center;'>💕💕Wedding Quiz💕💕</h1>",
        unsafe_allow_html=True,
    )

    def change_page():
        st.session_state.answers.append(st.session_state.answer0)
        st.session_state.page_id = "page1"

    with st.form("f0"):
        st.radio("テーブル番号を選んでね", ["A", "B", "C", "D", "E", "F", "G"], key="answer0")
        st.form_submit_button("スタート！", on_click=change_page)


# 問題１
def page1():
    st.markdown(
        "<h1 style='text-align: center;'>第１問</h1>",
        unsafe_allow_html=True,
    )

    def change_page():
        st.session_state.answers.append(st.session_state.answer1)
        st.session_state.page_id = "page2"

    with st.form("f1"):
        st.radio("新郎の名前は？", ["たろう", "じろう", "さぶろう"], key="answer1")
        st.form_submit_button("回答", on_click=change_page)


# 問題２
def page2():
    st.markdown(
        "<h1 style='text-align: center;'>第２問</h1>",
        unsafe_allow_html=True,
    )

    def change_page():
        st.session_state.answers.append(st.session_state.answer2)
        st.session_state.page_id = "page3"

    with st.form("f2"):
        st.radio("新婦の名前は？", ["はなこ", "ゆうこ", "うめこ"], key="answer2")
        st.form_submit_button("回答", on_click=change_page)


# 問題３
def page3():
    st.markdown(
        "<h1 style='text-align: center;'>第３問</h1>",
        unsafe_allow_html=True,
    )

    def change_page():
        st.session_state.answers.append(st.session_state.answer3)
        st.session_state.page_id = "page4"

    with st.form("f3"):
        st.radio("新郎の年齢は？", ["20歳", "30歳", "40歳", "50歳"], key="answer3")
        st.form_submit_button("回答", on_click=change_page)


# 問題４
def page4():
    st.markdown(
        "<h1 style='text-align: center;'>第４問</h1>",
        unsafe_allow_html=True,
    )

    def change_page():
        st.session_state.answers.append(st.session_state.answer4)
        st.session_state.page_id = "page5"

    with st.form("f4"):
        st.radio("新婦の年齢は？", ["20歳", "30歳", "40歳", "50歳"], key="answer4")
        st.form_submit_button("回答", on_click=change_page)


# 問題５
def page5():
    st.markdown(
        "<h1 style='text-align: center;'>第５問</h1>",
        unsafe_allow_html=True,
    )

    def change_page():
        st.session_state.answers.append(st.session_state.answer5)
        st.session_state.page_id = "page_end"

    with st.form("f5"):
        st.radio("２人の趣味は？", ["野球", "サッカー", "ピアノ", "卓球", "水泳"], key="answer5")
        st.form_submit_button("回答", on_click=change_page)


# 最終ページ
def page_end():

    # 回答内容をサーバに送ったりなんなり

    st.markdown(
        "<h1 style='text-align: center;'>回答ありがとう🎉</h1>",
        unsafe_allow_html=True,
    )
    st.markdown("---")
    st.markdown(
        "<h2 style='text-align: center;'>あなたの回答</h2>",
        unsafe_allow_html=True,
    )
    st.markdown(
        f"<div style='text-align: center;'>テーブル：{st.session_state.answers[0]}</div>",
        unsafe_allow_html=True,
    )
    for num, value in enumerate(st.session_state.answers, 0):
        if num != 0:
            st.markdown(
                f"<div style='text-align: center;'>第{num}問：{value}</div>",
                unsafe_allow_html=True,
            )
    ## バルーンを飛ばす
    st.balloons()


# ページ遷移のための判定
if st.session_state.page_id == "main":
    main()

if st.session_state.page_id == "page1":
    page1()

if st.session_state.page_id == "page2":
    page2()

if st.session_state.page_id == "page3":
    page3()

if st.session_state.page_id == "page4":
    page4()

if st.session_state.page_id == "page5":
    page5()

if st.session_state.page_id == "page_end":
    page_end()

