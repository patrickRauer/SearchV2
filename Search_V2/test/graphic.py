
from pycallgraph import PyCallGraph
from pycallgraph import Config
from pycallgraph import GlobbingFilter
from pycallgraph.output import GraphvizOutput
import test

if __name__ == '__main__':
    config = Config()
    config.trace_filter = GlobbingFilter(include=[
        'pycallgraph.*',
        'Search_V2.*',
        '*test*',
    ])
    graphviz = GraphvizOutput()
    graphviz.output_file = 'basic2.png'

    with PyCallGraph(output=graphviz, config=config):
        test.test_data_system()
