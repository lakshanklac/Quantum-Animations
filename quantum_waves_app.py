import streamlit as st

# Define two functions for the two different codes to be run
def run_first_code():
    st.write("Running the first code...")
    # Your first block of code logic here
    st.write("First code executed.")

def run_second_code():
    st.write("Running the second code...")
    # Your second block of code logic here
    st.write("Second code executed.")

# Create a title for the app
st.title("Webpage with Two Tiles")

# Create two buttons representing the two tiles
col1, col2 = st.columns(2)

with col1:
    if st.button("Run First Code"):
        run_first_code()

with col2:
    if st.button("Run Second Code"):
        run_second_code()
