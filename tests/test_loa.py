from desi_retriever.plotter import plot as desi_plot
from desi_retriever.loa import get_specs, get_rvspec_models


def test_loa_retrieval():
    SP = get_specs(survey='sv1',
                   program='dark',
                   hpx=17683,
                   targetid=39627652591521181)[0]
    assert SP is not None
    SP = get_specs(survey='sv1',
                   program='dark',
                   hpx=17683,
                   targetid=39627652591521181,
                   fibermap=True,
                   mask=True,
                   ivar=True)[0]
    assert SP is not None
    SPM = get_rvspec_models(survey='sv1',
                            program='dark',
                            hpx=17683,
                            targetid=39627652591521181)[0]
    assert SPM is not None
    my_sourceid = 1384229274732017408
    SP = get_specs(gaia_edr3_source_id=my_sourceid)[0]
    SPM = get_rvspec_models(gaia_edr3_source_id=my_sourceid)[0]
    desi_plot(SP, model=SPM)
    assert SP is not None
    assert SPM is not None
