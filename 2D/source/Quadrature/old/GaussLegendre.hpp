#include <vector>

#ifndef _GAUSSLEGENDRE_HPP
#define _GAUSSLEGENDRE_HPP

namespace Quadrature
{
    class GaussLegendre1D
    {
    public:
        GaussLegendre1D(unsigned doe);

        unsigned n_points() const;
        double weight(unsigned i) const;
        double point(unsigned i) const;

    private:
        const unsigned _doe;
        const unsigned _npts = _doe % 2 == 0 ? (_doe + 2) / 2 : (_doe + 1) / 2; // std::ceil((_doe + 1) / 2.0)

        std::vector<double> _points;
        std::vector<double> _weights;

        void sub_rule_01();
        void sub_rule_02();
        void sub_rule_03();
        void sub_rule_04();
        void sub_rule_05();
        void sub_rule_06();
        void sub_rule_07();
        void sub_rule_08();
        void sub_rule_09();
        void sub_rule_10();
        void sub_rule_11();
        void sub_rule_12();
        void sub_rule_13();
        void sub_rule_14();
        void sub_rule_15();
        void sub_rule_16();
        void sub_rule_17();
        void sub_rule_18();
        void sub_rule_19();
        void sub_rule_20();
        void sub_rule_21();
        void sub_rule_22();
        void sub_rule_23();
        void sub_rule_24();
        void sub_rule_25();
        void sub_rule_26();
        void sub_rule_27();
        void sub_rule_28();
        void sub_rule_29();
        void sub_rule_30();
        void sub_rule_31();
    };
}

#endif