import pandas as pd
from pathlib import Path

datapath = Path(input("Paste path to datafile"))
with open(datapath, 'r') as file:
    try:
        df = pd.read_csv(file, delimiter = '\t')
        try: 
            df = pd.read_csv(file, delimiter = ';')
        finally:
            print(df)