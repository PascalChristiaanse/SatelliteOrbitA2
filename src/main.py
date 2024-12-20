import numpy as np

from PartA import PartA
from PartB import PartB
from PartC import PartC


def main():
    # How large are the following effects on the pseudoranges measured by the GPS receiver?
    # 1) GPS clock offsets (4 points)
    # 2) Light time effect (10 points)
    # 3) Relativistic effect caused by the eccentricity of the GPS orbits (6 points)

    # 1.
    part_a = PartA()
    print("Q.A) Pseudo ranges")
    print("1)", part_a.Q1(), "m")
    print("2)", part_a.Q2(), "m")
    print("3)", np.abs(part_a.Q3()), "m")

    print("\nQ.B) PRN_ID File")
    part_b = PartB()
    print("1)", part_b.Q1())

    print("\nQ.C) Covariance matrix")
    part_c = PartC()
    print("1)", part_c.Q1())


if __name__ == "__main__":
    main()
