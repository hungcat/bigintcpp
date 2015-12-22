#include "bigint.h"

int main() {
    const uint_long len = 30;

    BigInt::initRandom();

    BigInt x = BigInt::randomLength(len);
    BigInt y = BigInt::randomBits(len / 2);

//    std::cout << x << std::endl;
//    std::cout << y << std::endl;
    BigInt a("4294324896967296");
    BigInt b("8862346823468364");
    BigInt c("5379398739324795");

    x = a * b;
    y = a * c;

    std::cout << a << std::endl;
    std::cout << x << std::endl;
    std::cout << y << std::endl;

    std::cout << "Gcd" << std::endl;
    std::cout << BigInt::basicGcd(x, y) << std::endl;
    std::cout << "BinaryGcd" << std::endl;
    std::cout << BigInt::binaryGcd(x, y) << std::endl;


/*
    BigInt z1 = BigInt::ZERO;
    BigInt z2(0);
    BigInt z3("0");

    std::cout << z1 << std::endl;
    std::cout << z2 << std::endl;
    std::cout << z3 << std::endl;
    z1.normalize();
    z2.normalize();
    z3.normalize();

    std::cout << z1 << std::endl;
    std::cout << z2 << std::endl;
    std::cout << z3 << std::endl;
    */
    /*
    BigInt hoge("12345678900987654360");
    BigInt hage("1234567850");
    // BigInt hoge("1000000000000000000000001");
    // BigInt hage("100000000000000000");
    std::cout << "hoge" << std::endl;
    std::cout << hoge.toStr() << std::endl;
    std::cout << "hage" << std::endl;
    std::cout << hage.toStr() << std::endl;
    */

    /*
    std::cout << hoge / hage << std::endl;
    std::cout << hoge % hage << std::endl;
    std::cout << hoge / hage.flip() << std::endl;
    std::cout << hoge % hage << std::endl;
    std::cout << hoge.flip() / hage << std::endl;
    std::cout << hoge % hage << std::endl;
    std::cout << hoge / hage.flip() << std::endl;
    std::cout << hoge % hage << std::endl;

    for (int i = 0; i < 20; i++) {
        std::cout << BigInt::factorial(i) << std::endl;
    }

    */




    /*
    // BigInt x = BigInt::randomLength(len);
    BigInt y = BigInt::randomBits(len);

    // std::cout << x << std::endl;
    std::cout << y << std::endl;

    BigInt z = y.sqrt();
    std::cout << z << std::endl;
    std::cout << y.div(z) << std::endl;

    BigInt z = BigInt::randomLength(len*2);
    std::cout << z << std::endl;

    std::cout << "divrem" << std::endl;
    BigInt quot, rem;
    BigInt::divrem(z, x, &quot, &rem);

    std::cout << quot << std::endl;
    std::cout << rem << std::endl;

    std::cout << "comp" << std::endl;
    std::cout << (quot * x + rem - z).isZero() << std::endl;
    */

    return 0;
}
