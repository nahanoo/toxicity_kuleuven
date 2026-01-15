import plotly.graph_objects as go
import pandas as pd
from os.path import join
from lmfit.model import load_modelresult


m = load_modelresult("models/ST_1.lmfit")
