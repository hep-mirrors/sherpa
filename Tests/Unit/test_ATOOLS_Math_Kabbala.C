#include <catch2/catch_all.hpp>
#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Org/Message.H"
#include <map>
#include <cmath>

using namespace ATOOLS;

TEST_CASE("KabbalaConstructors", "[ATOOLS::Math::Kabbala][ATOOLS][Math]"){
    /*INIT*/
    std::map<std::string, double_t>* mp = new std::map<std::string, double_t>();
    mp->insert(std::make_pair("a", 2.));
    mp->insert(std::make_pair("b", 4.));
    mp->insert(std::make_pair("c", 0.125));
    /*TESTS*/
    SECTION("Default"){
        Kabbala k_0 = Kabbala();
        CHECK(k_0.Value() == Complex(0.,0.));
        CHECK(k_0.String() == "0");
        k_0.Update(mp);
        CHECK(k_0.Value() == Complex(0.,0.));
        CHECK(k_0.String() == "0");
    }
    SECTION("ComplexOnly"){
        Kabbala k_1 = Kabbala(1.0);
        CHECK(k_1.Value() == Complex(1.,0.));
        CHECK(k_1.String() == ToString(1.));
        k_1.Update(mp);
        CHECK(k_1.Value() == Complex(1.,0.));
        CHECK(k_1.String() == ToString(1.));
    }
    SECTION("Legacy"){
        Kabbala k_a = Kabbala("a", Complex(2., -1.));
        CHECK(k_a.Value() == Complex(2., -1.));
        CHECK(k_a.String() == "a");
        k_a.Update(mp);
        CHECK(k_a.Value() == Complex(2., 0.));
        CHECK(k_a.String() == "a");
    }
    SECTION("SpecifyAll"){
        Kabbala::Func fb = [](Kabbala::Function_Argument map){return Complex(0,4);};
        Kabbala k_b = Kabbala("b", 2.0, fb);
        CHECK(k_b.Value() == Complex(2., 0.));
        CHECK(k_b.String() == "b");
        k_b.Update(mp);
        CHECK(k_b.Value() == Complex(0., 4.));
        CHECK(k_b.String() == "b");
        Kabbala::Func fdb = [](Kabbala::Function_Argument map){if (map->count("b")==0) return C_ZERO; return Complex(map->at("b"));};
        Kabbala k_db = Kabbala("d", -1., fdb);
        CHECK(k_db.Value() == Complex(-1., 0.));
        CHECK(k_db.String() == "d");
        k_db.Update(mp);
        CHECK(k_db.Value() == Complex(4., 0.));
        CHECK(k_db.String() == "d");
    }
    SECTION("ValueFromMap"){
        Kabbala::Func fc = [](Kabbala::Function_Argument map){if (map->count("c")==0) return C_ZERO; return Complex(map->at("c"));};
        Kabbala k_c = Kabbala("c", fc, mp);
        CHECK(k_c.Value() == Complex(0.125, 0.));
        CHECK(k_c.String() == "c");
        k_c.Update(mp);
        CHECK(k_c.Value() == Complex(0.125, 0.));
        CHECK(k_c.String() == "c");
    }
    /*CLEANUP*/
    delete mp;
}

TEST_CASE("KabbalaOperations", "[ATOOLS::Math::Kabbala][ATOOLS][Math]"){
    /*INIT*/
    std::map<std::string, double_t>* mp = new std::map<std::string, double_t>();
    mp->insert(std::make_pair("a", 2.));
    mp->insert(std::make_pair("b", 4.));
    mp->insert(std::make_pair("c", 0.125));
    Kabbala k1, k2, k3, k;
    /*TESTS*/
    SECTION("Negation"){
        k1 = Kabbala("a", Complex(1., 0.));
        k2 = -k1;
        CHECK(k2.Value() == Complex(-1., 0.));
        CHECK(k2.String() == "-(a)");
        k2.Update(mp);
        CHECK(k2.Value() == Complex(-2., 0.));
    }
        
    SECTION("Addition and Subtraction") {
        k1 = Kabbala("a", Complex(1., 0.));
        k2 = Kabbala("b", Complex(0., -1.));
        k3 = k1 + k2;
        CHECK(k3.Value() == Complex(1., -1.));
        CHECK(k3.String() == "a+b");
        k3.Update(mp);
        CHECK(k3.Value() == Complex(6., 0.));

        k = k1 - k2;
        CHECK(k.Value() == Complex(1., 1.));
        CHECK(k.String() == "a-(b)");
        k.Update(mp);
        CHECK(k.Value() == Complex(-2., 0.));
    }

    SECTION("In-place Addition and Subtraction") {
        k1 = Kabbala("a", Complex(1., 0.));
        k2 = Kabbala("b", Complex(0., -1.));
        k1 += k2;
        CHECK(k1.Value() == Complex(1., -1.));
        CHECK(k1.String() == "a+b");
        k1.Update(mp);
        CHECK(k1.Value() == Complex(6., 0.));

        k1 -= k2;
        CHECK(k1.Value() == Complex(6., 1.));
        CHECK(k1.String() == "a+b-(b)");
        k1.Update(mp);
        CHECK(k1.Value() == Complex(2., 0.));
    }

    SECTION("Multiplication") {
        k1 = Kabbala("a", Complex(2., 0.));
        k2 = Kabbala("b", Complex(0., -1.));
        k3 = k1 * k2;
        CHECK(k3.Value() == Complex(0., -2.));
        CHECK(k3.String() == "(a)*(b)");
        k3.Update(mp);
        CHECK(k3.Value() == Complex(8., 0.));
    }

    SECTION("Division") {
        k1 = Kabbala("a", Complex(4., 0.));
        k2 = Kabbala("b", Complex(2., 0.));
        k3 = k1 / k2;
        CHECK(k3.Value() == Complex(2., 0.));
        CHECK(k3.String() == "(a)/(b)");
        k3.Update(mp);
        CHECK(k3.Value() == Complex(0.5, 0.));
    }

    SECTION("Exponentiation") {
        k1 = Kabbala("a", Complex(0.5, 0.));
        k2 = pow(k1, Complex(2., 0.));
        CHECK(k2.Value() == Complex(0.25, 0.));
        CHECK(k2.String() == "(a)^(2)");
        k2.Update(mp);
        CHECK(k2.Value() == Complex(4., 0.));
    }

    SECTION("Square Root") {
        k1 = Kabbala("b", Complex(9., 0.));
        k2 = sqrt(k1);
        CHECK(k2.Value() == Complex(3., 0.));
        CHECK(k2.String() == "sqrt(b)");
        k2.Update(mp);
        CHECK(k2.Value() == Complex(2., 0.));

        // Check square root of negative
        k1 = Kabbala("b", Complex(-4., 0.));
        k2 = sqrt(k1);
        CHECK(k2.Value() == Complex(0., 2.));
        CHECK(k2.String() == "sqrt(b)");
    }

    SECTION("Trigonometric Functions") {
        k1 = Kabbala("a", Complex(0., 0.5));
        
        // Sin
        k2 = sin(k1);
        CHECK(k2.Value() == std::sin(Complex(0., 0.5)));
        CHECK(k2.String() == "sin(a)");
        k2.Update(mp);
        CHECK(k2.Value() ==  std::sin(Complex(2., 0.)));

        // Cos
        k2 = cos(k1);
        CHECK(k2.Value() == std::cos(Complex(0., 0.5)));
        CHECK(k2.String() == "cos(a)");
        k2.Update(mp);
        CHECK(k2.Value() == std::cos(Complex(2., 0.)));

        // Tan
        k2 = tan(k1);
        CHECK(k2.Value() == std::tan(Complex(0., 0.5)));
        CHECK(k2.String() == "tan(a)");
        k2.Update(mp);
        CHECK(k2.Value() == std::tan(Complex(2., 0.)));
    }

    SECTION("Complex Conjugate") {
        k1 = Kabbala("a", Complex(2., 3.));
        k2 = complexconjugate(k1);
        CHECK(k2.Value() == Complex(2., -3.));
        CHECK(k2.String() == "(a)*");
        k2.Update(mp);
        CHECK(k2.Value() == Complex(2., 0.));
    }
    /*CLEANUP*/
    delete mp;
}
