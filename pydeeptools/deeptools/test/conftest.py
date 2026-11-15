import matplotlib
import matplotlib.pyplot as plt
import pytest


@pytest.fixture(autouse=True)
def _isolate_matplotlib_rcparams():
    with matplotlib.rc_context():
        yield
    plt.close("all")
