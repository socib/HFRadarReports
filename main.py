from report_utils import *


def main():
    print get_thredds_opendap_link('hf_radar', 'hf_radar_ibiza-scb_codarssproc001', 1, 'dep0001', 2016, 4)


if __name__ == "__main__":
    main()
