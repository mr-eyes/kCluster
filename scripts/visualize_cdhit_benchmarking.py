import re
import sys
import os
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.offline

output_html = "visualization.html"

if len(sys.argv) > 1:
    output_html = sys.argv[1]

thresholds = []
types = {}

#kmers_count = os.popen('grep "" kmers_clustering/full_human_transcriptome_results/*/*summary*').read()

#for line in kmers_count.split("\n"):
for line in sys.stdin:
    if len(line) < 5:
        continue

    _threshold = re.findall(r"(\d+)%", line)[0]
    _value = re.findall(r":(\d+)", line)[0]
    _type = re.findall(r"\[_(\w+)\]", line)[0]

    if _threshold not in thresholds:
        thresholds.append(_threshold)
        
    if _type not in types:
        types[_type] = [_value]
    else:
        types[_type].append(_value)


complete_clean = go.Scatter(
    x = thresholds,
    y = types["complete_clean"],
    mode = 'lines+markers',
    name = 'Complete Clean'
)

complete_mixed = go.Scatter(
    x = thresholds,
    y = types["complete_mixed"],
    mode = 'lines+markers',
    name = 'Complete Mixed'
)

incomplete_clean = go.Scatter(
    x = thresholds,
     y = types["incomplete_clean"],
    mode = 'lines+markers',
    name = 'Incomplete Clean'
)

incomplete_mixed = go.Scatter(
    x = thresholds,
    y = types["incomplete_mixed"],
    mode = 'lines+markers',
    name = 'Incomplete Mixed'
)

data = [complete_mixed, complete_clean, incomplete_mixed, incomplete_clean]

layout = dict(title = '<b>kSpider vs CD-HIT Clustering Assessment</b>',
              xaxis = dict(title = '<b>Threshold</b>', type='category'),
              yaxis = dict(title = '<b>Number Of Clusters</b>', type="log"),
              )

fig = dict(data=data, layout=layout)

plotly.offline.plot(fig, filename = output_html, auto_open = False)
