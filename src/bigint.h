#ifndef HUNGCAT_BIGINT_CPP
#define HUNGCAT_BIGINT_CPP

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <random>
#include <sstream>

using int_long = std::int_fast64_t;
using uint_long = std::uint_fast64_t;
using int_rep = std::int32_t;
using uint_rep = std::uint32_t;

class BigInt;
class MontgomerySystem;

/*
 *
 *
 * Definition   {{{
 *
 *
 */

/*
 *
 * BigInt {{{
 *
 */

class BigInt {

    /*
     *
     * Privates
     *
     */

    /*
     * Class constants
     */

    static const uint_long BIT_SIZE;    // uint_repのビット長
    static const uint_long LOWER_MASK;  // BIT_SIZE bit lower mask
    static const uint_long UPPER_MASK;  // BIT_SIZE bit upper mask

    static const uint_long DEC_SIZE;  // uint_repの10進長(切り捨て)
    static const BigInt BASE10;  // 2^BIT_SIZE未満の最大の10のべき乗

    static const std::string NUMBER_LIST;  // 数値表現用リスト

    /*
     * Class variables
     */

    static std::mt19937 mt32;     // 32bit mt
    static std::mt19937_64 mt64;  // 64bit mt
    static MontgomerySystem _M;

    /*
     * Member variables
     */

    bool _NaN;                   // true: not a number
    bool _sign;                  // true: negative
    std::vector<uint_rep> _rep;  // internal representation

    /*
     * Member methods
     */

    BigInt& pushUpper(uint_long n);     // 内部配列の上位桁に追加
    BigInt& pushLower(uint_long n);     // 内部配列の下位桁に追加
    BigInt& normalize();                // 正規化
    BigInt& NaN();                      // NaN化
    BigInt& Zero();                     // ゼロ化
    BigInt& shift(int_long num);        // 内部配列単位でシフト
    BigInt& ushift(uint_long num = 1);  // 左シフト(上位桁方向)
    BigInt& lshift(uint_long num = 1);  // 右シフト(下位桁方向)

public:
    /*
     *
     * Publics
     *
     */

    /*
     * Constructor
     */

    // 必ずこいつが呼ばれる
    BigInt(uint_long n) : _NaN(false), _sign(false), _rep() { pushLower(n); }
    BigInt() : BigInt(0UL) {}
    BigInt(int_long n) : BigInt((uint_long)(n < 0 ? -n : n)) {
        _sign = (n < 0);
    }
    BigInt(uint_rep n) : BigInt((uint_long)n) {}
    BigInt(int_rep n) : BigInt((int_long)n) {}
    BigInt(std::string str, uint_long radix = 10) : BigInt() {
        BigInt::parse(*this, str, radix);
    }

    /*
     * CONSTANTS
     */

    static const BigInt ZERO;
    static const BigInt ONE;
    static const BigInt TWO;

    /*
     * General
     */

    BigInt& setSize(uint_long s);  // 内部配列の長さをセット
    BigInt& flip();                // 正負を切り替える

    static uint_long bitCount(uint_long num);  // 立っているビットの数
    template <typename T>
    static uint_long numOfLeadingZeros(T num);  // MSBからの0の数
    template <typename T>
    static uint_long numOfTrailingZeros(T num);  // LSBからの0の数
    template <typename T>
    static uint_long bitLength(T num);  // bit長

    uint_long bitCount() const;            // 立っているビットの数
    uint_long numOfLeadingZeros() const;   // MSBからの0の数
    uint_long numOfTrailingZeros() const;  // LSBからの0の数
    uint_long bitLength() const;           // bit長

    BigInt& set(uint_long i);
    BigInt& unset(uint_long i);

    /*
     * Bit Arithmetic
     */

    BigInt& operator&=(const BigInt& o);
    BigInt& operator|=(const BigInt& o);
    BigInt& operator^=(const BigInt& o);
    BigInt& operator<<=(uint_long num);  // 左シフト(bit単位)
    BigInt& operator>>=(uint_long num);  // 右シフト(bit単位)

    BigInt operator&(const BigInt& o) const;
    BigInt operator|(const BigInt& o) const;
    BigInt operator^(const BigInt& o) const;
    BigInt operator<<(uint_long num) const;  // 左シフト(bit単位)
    BigInt operator>>(uint_long num) const;  // 右シフト(bit単位)

    /*
     * Comparison
     */

    int_long comp(const BigInt& o) const;     // 比較
    int_long compAbs(const BigInt& o) const;  // 絶対値比較

    bool operator<(const BigInt& o) const { return comp(o) < 0; }
    bool operator>(const BigInt& o) const { return comp(o) > 0; }
    bool operator==(const BigInt& o) const { return comp(o) == 0; }
    bool operator<=(const BigInt& o) const { return !operator>(o); }
    bool operator>=(const BigInt& o) const { return !operator<(o); }
    bool operator!=(const BigInt& o) const { return !operator==(o); }

    bool isNaN() const { return _NaN; }
    bool isZero() const { return !isNaN() && operator==(ZERO); }
    bool isPositive() const { return !isNaN() && operator>(ZERO); }
    bool isNegative() const { return !isNaN() && operator<(ZERO); }
    bool isOdd() const { return !isNaN() && _rep.size() > 0 && (_rep[0] & 1); }
    bool isEven() const {
        return !isNaN() && _rep.size() > 0 && !(_rep[0] & 1);
    }

    /*
     * Basic Arithmetic
     */

    BigInt& operator+=(const BigInt& o);
    BigInt& operator-=(const BigInt& o);
    BigInt& operator*=(const BigInt& o);
    BigInt& operator/=(const BigInt& o);
    BigInt& operator%=(const BigInt& o);
    BigInt operator+(const BigInt& o) const { return add(o); }
    BigInt operator-(const BigInt& o) const { return sub(o); }
    BigInt operator*(const BigInt& o) const { return mul(o); }
    BigInt operator/(const BigInt& o) const { return div(o); }
    BigInt operator%(const BigInt& o) const { return rem(o); }

    BigInt& operator*=(uint_rep o);
    BigInt& operator*=(int_rep o);
    BigInt& operator*=(uint_long o);
    BigInt& operator*=(int_long o);
    BigInt operator*(uint_rep o) const { return BigInt(*this).operator*=(o); }
    BigInt operator*(int_rep o) const { return BigInt(*this).operator*=(o); }
    BigInt operator*(uint_long o) const { return BigInt(*this).operator*=(o); }
    BigInt operator*(int_long o) const { return BigInt(*this).operator*=(o); }

    BigInt add(const BigInt& o) const;
    BigInt sub(const BigInt& o) const;

    BigInt mul(const BigInt& o) const;
    BigInt mulTrivial(const BigInt& o) const;

    BigInt div(const BigInt& o, BigInt* rem = NULL) const;
    BigInt rem(const BigInt& o, BigInt* div = NULL) const;
    static bool divrem(const BigInt& numer, const BigInt& denom, BigInt* q,
                       BigInt* r);
    static bool divremTrivial(const BigInt& numer, const BigInt& denom,
                              BigInt* q, BigInt* r);

    /*
     * Other computations
     */

    static BigInt mulMod(const BigInt& a, const BigInt& b, const BigInt& n);
    static BigInt mulMod(const BigInt& a, const BigInt& b,
                         const MontgomerySystem& M);
    static BigInt powMod(const BigInt& a, uint_long e, const BigInt& n);
    static BigInt powMod(const BigInt& a, uint_long e,
                         const MontgomerySystem& M);
    BigInt mulMod(const BigInt& o, const BigInt& n);
    BigInt mulMod(const BigInt& o, const MontgomerySystem& M);
    BigInt powMod(uint_long e, const BigInt& n);
    BigInt powMod(uint_long e, const MontgomerySystem& M);

    BigInt sqr() const { return mul(*this); }
    BigInt pow(uint_long exp) const;
    BigInt pow(const BigInt& exp) const;

    static uint_rep sqrt(uint_long x);
    BigInt sqrt() const;
    BigInt sqrtTrivial() const;
    BigInt sqrtNewton() const;

    BigInt factorial() const;
    static BigInt factorial(uint_long n) { return BigInt(n).factorial(); }

    static BigInt gcd(const BigInt& a, const BigInt& b);
    static BigInt basicGcd(const BigInt& a, const BigInt& b);
    static BigInt binaryGcd(const BigInt& a, const BigInt& b);

    /*
     * Parser
     * Radix Conversion
     */

    static bool parse(BigInt& bint, std::string str, uint_long radix = 10);
    static BigInt parse(const std::string& str, uint_long radix = 10);
    static std::string toStr(const BigInt& bint, uint_long radix = 10);
    static std::string uintToStr(uint_long num, uint_long radix = 10);
    std::string toStr(uint_long radix = 10) const;
    std::string toStrAsVector() const;
    friend std::ostream& operator<<(std::ostream& os, const BigInt& bint);

    /*
     * Randoms
     */

    // cryptおぷしょんはみじっそう
    static uint_long initRandom(uint_long s = 0);
    static BigInt randomBound(uint_long n, bool crypt = false);
    static BigInt randomBound(const BigInt& n, bool crypt = false);
    static BigInt randomBits(uint_long n, bool crypt = false);
    static BigInt randomLength(uint_long n, bool crypt = false);
};
/*
 *
 * }}}
 *
 */

/*
 *
 * MontgomerySystem {{{
 *
 */

class MontgomerySystem {
    BigInt N;        // modulus(odd)
    BigInt Np;       // N * Np == -1 mod R
    BigInt R2;       // R*R mod N
    BigInt R_MASK;   // R - 1
    BigInt NR;       // N*R: limit of arg to reduce
    uint_long lenR;  // bitlength of R
                     // BigInt R;   // gcd(R, N) == 1, R = 2^hoge

    bool isValid(const BigInt& x) const { return (!x.isNegative() && x < NR); }

public:
    MontgomerySystem() : N(), Np(), R2(), R_MASK(), NR(), lenR() {}
    MontgomerySystem(const BigInt& n) : MontgomerySystem() { setParams(n); }

    bool isSet() const { return !N.isZero(); }
    const BigInt& getModulus() const { return N; }
    const BigInt& getR2() const { return R2; }

    bool setParams(const BigInt& n) {
        if (n.isEven() || n.isNegative()) return false;

        N = n;
        uint_long nlen = N.bitLength();
        BigInt R = BigInt::ONE;
        R <<= nlen;
        R2 = (R * R) % N;
        R_MASK = R - 1;
        lenR = nlen;
        NR = N << lenR;

        // N*Np ≡ -1 mod R
        Np = 0; /* 求めるN' */
        BigInt t = 0;
        for (uint_long i = 0; i < nlen; ++i) {
            if (t.isEven()) {
                /* ゼロになっているビットがあったら、N'のその部分を1にする（NはRと互いに素なので必ず奇数）*/
                t += N; /* 掛け算だが、二進数一桁の掛け算なので実質は足し算 */
                Np.set(i); /* N'のその部分を1にする */
            }
            t >>= 1; /* 必ず端数が出るが切り捨てる */
        }

        //      return true;
        return (((N * Np) + BigInt::ONE) % R).isZero();
    }

    BigInt& _reduce(BigInt& T) const {
        if (!isValid(T)) {
            T %= N;
            if (T.isNegative()) T += N;
        }
        BigInt t = T;
        T.operator*=(Np)
            .operator&=(R_MASK)
            .operator*=(N)
            .operator+=(t)
            .operator>>=(lenR);
        // T = (((T * Np) & R_MASK) * N + T) >> lenR;
        return (T >= N) ? T.operator-=(N) : T;
    }
    BigInt reduce(const BigInt& T) const {
        BigInt t(T);
        _reduce(t);
        return t;
    }
    BigInt& _toME(BigInt& T) const { return _reduce(T.operator*=(R2)); }
    BigInt toME(const BigInt& T) const {
        BigInt t(T);
        _toME(t);
        return t;
    }
};

/*
 *
 * }}}
 *
 */

/*
 *
 *
 * }}}
 *
 *
 */

/*
 *
 * Static members {{{
 *
 */

// private constatns
const uint_long BigInt::BIT_SIZE =
    std::numeric_limits<uint_rep>::digits;          // uint_repのビット長
const uint_long BigInt::LOWER_MASK = ~(uint_rep)0;  // BIT_SIZE bit lower mask
const uint_long BigInt::UPPER_MASK = ~LOWER_MASK;   // BIT_SIZE bit upper mask

const uint_long BigInt::DEC_SIZE =
    BIT_SIZE * std::log10(2);  // uint_repの10進長
const BigInt BigInt::BASE10 =
    BigInt(10).pow(BigInt::DEC_SIZE);  // 2^BIT_SIZE未満の最大の10のべき乗

const std::string BigInt::NUMBER_LIST(
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");  // 数値表現用リスト

// public constants
const BigInt BigInt::ZERO(0);
const BigInt BigInt::ONE(1);
const BigInt BigInt::TWO(2);

std::mt19937 BigInt::mt32;     // 32bit mt
std::mt19937_64 BigInt::mt64;  // 64bit mt
MontgomerySystem BigInt::_M;   // for Montgomery Arithmetic

/*
 *
 * }}}
 *
 */

/*
 *
 * Private member methods {{{
 *
 */

/*
 * General {{{
 */

// 上位桁に追加(BIT_SIZE進数)
// (0の場合は追加しない)
inline BigInt& BigInt::pushUpper(uint_long n) {
    if (isNaN()) return *this;

    if (n) {
        _rep.push_back(n);  //  & LOWER_MASK
        if (n & UPPER_MASK) _rep.push_back(n >> BIT_SIZE);
    }

    return normalize();
}
// 下位桁に追加(BIT_SIZE進数)
inline BigInt& BigInt::pushLower(uint_long n) {
    if (isNaN()) return *this;

    auto b = begin(_rep);
    if (n & UPPER_MASK) b = _rep.insert(b, n >> BIT_SIZE);
    _rep.insert(b, n);  //  & LOWER_MASK

    return normalize();
}
// 正規化
inline BigInt& BigInt::normalize() {
    if (isNaN()) return NaN();

    uint_long max = _rep.size();

    if (max == 0) return Zero();

    while (--max && _rep[max] == 0)
        ;
    setSize(max + 1);
    if (isZero()) _sign = false;
    return *this;
}
// NaN化
inline BigInt& BigInt::NaN() {
    _NaN = true;
    _sign = false;
    _rep.clear();
    return *this;
}
// ゼロ化
inline BigInt& BigInt::Zero() {
    _NaN = false;
    _sign = false;
    _rep.assign(1, 0);
    return *this;
}
// _rep shift(BIT_SIZE進数)
inline BigInt& BigInt::shift(int_long num) {
    if (num > 0) ushift(num);
    if (num < 0) lshift(-num);
    return *this;
}
// _rep upper shift(BIT_SIZE進数)
inline BigInt& BigInt::ushift(uint_long num) {
    if (isNaN() || isZero()) return *this;

    _rep.insert(begin(_rep), num, 0);
    return *this;
}
// _rep lower shift(BIT_SIZE進数)
inline BigInt& BigInt::lshift(uint_long num) {
    if (isNaN() || isZero()) return *this;

    auto b = begin(_rep);
    _rep.erase(b, b + num);
    return *this;
}

/*
 * }}}
 */

/*
 *
 * }}}
 *
 */

/*
 *
 * Public member methods {{{
 *
 */

/*
 * General {{{
 */

// 桁数指定(BIT_SIZE進数)
inline BigInt& BigInt::setSize(uint_long s) {
    if (isNaN()) return *this;

    _rep.resize(s, 0);
    return *this;
}
// sign反転
inline BigInt& BigInt::flip() {
    if (isNaN()) return *this;

    _sign = !_sign;
    return *this;
}

// 64bit counter
inline uint_long BigInt::bitCount(uint_long num) {
    num -= num >> 1 & 0x5555555555555555;
    num = (num & 0x3333333333333333) + (num >> 2 & 0x3333333333333333);
    num = (num + (num >> 4)) & 0x0F0F0F0F0F0F0F0F;
    num += num >> 8;
    num += num >> 16;
    return (num + (num >> 32)) & 0xFF;
}

// numのMSBから続く0の長さを返す
// bitlength := bitsize(T) - NLZ
template <typename T>
inline uint_long BigInt::numOfLeadingZeros(T num) {
    const uint_long bitsize = std::numeric_limits<T>::digits;

    if (num == 0) return bitsize;

    uint_long c = 0, t;

    t = num & 0xFFFFFFFF00000000;
    if (t) {
        num = t;
        c |= 0x20;
    }
    t = num & 0xFFFF0000FFFF0000;
    if (t) {
        num = t;
        c |= 0x10;
    }
    t = num & 0xFF00FF00FF00FF00;
    if (t) {
        num = t;
        c |= 0x08;
    }
    t = num & 0xF0F0F0F0F0F0F0F0;
    if (t) {
        num = t;
        c |= 0x04;
    }
    t = num & 0xCCCCCCCCCCCCCCCC;
    if (t) {
        num = t;
        c |= 0x02;
    }
    if (num & 0xAAAAAAAAAAAAAAAA) {
        c |= 0x01;
    }
    return c ^ (bitsize - 1);
}

// numのLSBから続く0の長さを返す
template <typename T>
inline uint_long BigInt::numOfTrailingZeros(T num) {
    if (num == 0) return std::numeric_limits<T>::digits;
    return BigInt::bitCount((~num) & (num - 1));
}
// numのビット長を返す(最大64bit)
// floor(log_2{num}) := bitlength - 1
// NLZ := 64 - bitlength
template <typename T>
inline uint_long BigInt::bitLength(T num) {
    return std::numeric_limits<T>::digits - BigInt::numOfLeadingZeros(num);
}
// numのMSBから続く0の長さを返す
// 配列が過剰に用意されている場合であっても、0埋めされているとしてカウントする
inline uint_long BigInt::numOfLeadingZeros() const {
    if (isNaN()) return -1;

    uint_long s = _rep.size(), i = s;
    if (s == 0) return 0;
    while (--i && _rep[i] == 0)
        ;
    return ((s - 1) - i) * BIT_SIZE + numOfLeadingZeros(_rep[i]);
}
// numのLSBから続く0の長さを返す
inline uint_long BigInt::numOfTrailingZeros() const {
    if (isNaN()) return -1;

    uint_long s = _rep.size(), i = 0;
    if (s == 0) return 0;
    for (i = 0; i < s - 1 && _rep[i] == 0; ++i)
        ;
    return i * BIT_SIZE + numOfTrailingZeros(_rep[i]);
}
// 自身のビット長を返す
inline uint_long BigInt::bitLength() const {
    if (isNaN()) return -1;

    auto s = _rep.size();
    if (s == 0) return 0;
    while (--s && _rep[s] == 0)
        ;
    return s * BIT_SIZE + bitLength(_rep[s]);
}

inline BigInt& BigInt::set(uint_long i) {
    if (isNaN()) return *this;

    uint_long s = _rep.size();
    uint_long q = i / BIT_SIZE, r = i % BIT_SIZE;

    if (q < s) _rep[q] &= 1UL << r;

    return *this;
}
inline BigInt& BigInt::unset(uint_long i) {
    if (isNaN()) return *this;

    uint_long s = _rep.size();
    uint_long q = i / BIT_SIZE, r = i % BIT_SIZE;

    if (q < s) _rep[q] &= ~(1UL << r);

    return *this;
}

/*
 * }}}
 */

/*
 * Bit Arithmetic {{{
 */

inline BigInt& BigInt::operator&=(const BigInt& o) {
    if (o.isNaN()) NaN();
    if (isNaN()) return *this;

    uint_long s = _rep.size(), os = o._rep.size();
    if (s < os) {
        setSize(os);
        s = os;
    }

    while (s--) _rep[s] &= o._rep[s];

    return *this;
}
inline BigInt& BigInt::operator|=(const BigInt& o) {
    if (o.isNaN()) NaN();
    if (isNaN()) return *this;

    uint_long s = _rep.size(), os = o._rep.size();
    if (s < os) {
        setSize(os);
        s = os;
    }

    while (s--) _rep[s] |= o._rep[s];

    return *this;
}
inline BigInt& BigInt::operator^=(const BigInt& o) {
    if (o.isNaN()) NaN();
    if (isNaN()) return *this;

    uint_long s = _rep.size(), os = o._rep.size();
    if (s < os) {
        setSize(os);
        s = os;
    }

    while (s--) _rep[s] ^= o._rep[s];

    return *this;
}

inline BigInt& BigInt::operator<<=(uint_long num) {
    if (isNaN()) return *this;

    uint_long r = num % BIT_SIZE, r_ = BIT_SIZE - r;
    uint_long s = _rep.size() - 1;

    // 型のbit長でmoduloをとられた分だけシフトする(32とかは注意)
    if (r) {
        if (bitLength(_rep[s]) > r_) pushUpper(_rep[s] >> r_);
        if (s) {
            while (s--) {
                // _rep[s + 1] = ((_rep[s + 1] << r) & LOWER_MASK) | _rep[s]
                _rep[s + 1] = _rep[s + 1] << r | _rep[s] >> r_;
            }
        }
        //	_rep[0] = (_rep[0] << r) & LOWER_MASK;
        _rep[0] <<= r;
    }

    return ushift(num / BIT_SIZE);
}
inline BigInt& BigInt::operator>>=(uint_long num) {
    if (isNaN()) return *this;

    lshift(num / BIT_SIZE);

    uint_long r = num % BIT_SIZE, r_ = BIT_SIZE - r;
    uint_long size = _rep.size() - 1;

    // 型のbit長でmoduloをとられた分だけシフトする(32とかは注意)
    if (r) {
        for (uint_long s = 0; s < size; s++) {
            _rep[s] = _rep[s + 1] << r_ | _rep[s] >> r;
        }
        _rep[size] >>= r;
    }

    return normalize();
}

inline BigInt BigInt::operator&(const BigInt& o) const {
    return BigInt(*this).operator&=(o);
}
inline BigInt BigInt::operator|(const BigInt& o) const {
    return BigInt(*this).operator|=(o);
}
inline BigInt BigInt::operator^(const BigInt& o) const {
    return BigInt(*this).operator^=(o);
}
inline BigInt BigInt::operator<<(uint_long num) const {
    return BigInt(*this).operator<<=(num);
}
inline BigInt BigInt::operator>>(uint_long num) const {
    return BigInt(*this).operator>>=(num);
}

/*
 * }}}
 */

/*
 * Comparison {{{
 */

// 比較
inline int_long BigInt::comp(const BigInt& o) const {
    int_long sign = _sign ? -1 : 1;
    return _sign != o._sign ? sign : sign * compAbs(o);
}
// 絶対値比較
inline int_long BigInt::compAbs(const BigInt& o) const {
    if (isNaN() || o.isNaN()) return isNaN() ? o.isNaN() ? 0 : -1 : 1;

    int_long diff = bitLength() - o.bitLength();
    if (diff) return diff;
    uint_long s = _rep.size(), os = o._rep.size();
    if (s * os == 0) return 0;
    while (s-- && !(diff = (int_long)_rep[s] - (int_long)o._rep[s]))
        ;
    return diff;
}

/*
 * }}}
 */

/*
 * Addition & Subtraction {{{
 */

// 加算
inline BigInt& BigInt::operator+=(const BigInt& o) {
    if (o.isNaN()) NaN();
    if (isNaN()) return *this;

    if (_sign != o._sign) {
        return flip().operator-=(o).flip();
    }

    uint_long s = _rep.size(), os = o._rep.size();
    if (s < os) {
        setSize(os);
        s = os;
    }

    uint_long carry = 0;
    for (uint_long i = 0; i < s; i++) {
        if (i < os) {
            carry += o._rep[i];
        } else if (carry == 0) {
            break;
        }
        carry += _rep[i];
        _rep[i] = carry;  // & LOWER_MASK;
        carry = (carry & UPPER_MASK) ? 1 : 0;
    }

    return pushUpper(carry);
}
// 加算
inline BigInt BigInt::add(const BigInt& o) const {
    return BigInt(*this).operator+=(o);
}

// 減算
inline BigInt& BigInt::operator-=(const BigInt& o) {
    if (o.isNaN()) NaN();
    if (isNaN()) return *this;

    if (_sign != o._sign) {
        return flip().operator+=(o).flip();
    }

    uint_long s = _rep.size(), os = o._rep.size();
    if (s < os) {
        setSize(os);
        s = os;
    }

    uint_long carry = 0;
    for (uint_long i = 0; i < s; i++) {
        if (i < os) {
            carry += o._rep[i];
        } else if (carry == 0) {
            break;
        }
        carry = (uint_long)_rep[i] - carry;
        _rep[i] = carry;  // & LOWER_MASK;
        carry = (carry & UPPER_MASK) ? 1 : 0;
    }

    if (carry) {
        for (auto&& r : _rep) {
            r = ~r;
        }
        operator+=(ONE);
        flip();
    }

    return normalize();
}
// 減算
inline BigInt BigInt::sub(const BigInt& o) const {
    if (this == &o) return ZERO;
    return BigInt(*this).operator-=(o);
}
/*
 * }}}
 */

/*
 * Multiplication {{{
 */

// 乗算
inline BigInt& BigInt::operator*=(const BigInt& o) { return (*this = mul(o)); }
// 乗算 あとでKaratubaも書く
inline BigInt BigInt::mul(const BigInt& o) const { return mulTrivial(o); }
inline BigInt& BigInt::operator*=(int_rep o) { return operator*=((int_long)o); }
inline BigInt& BigInt::operator*=(int_long o) {
    if (o < 0) {
        o = -o;
        flip();
    }
    return operator*=((uint_long)o);
}
inline BigInt& BigInt::operator*=(uint_rep o) {
    return operator*=((uint_long)o);
}
inline BigInt& BigInt::operator*=(uint_long o) {
    if (isNaN()) return NaN();
    if (isZero() || o == 0) return Zero();

    if (o & UPPER_MASK) return operator*=(BigInt(o));

    uint_long carry = 0;
    for (auto&& v : _rep) {
        carry += v * o;
        v = carry;  // & LOWER_MASK;
        carry >>= BIT_SIZE;
    }
    pushUpper(carry);

    return *this;
}
// trivial multiplication
inline BigInt BigInt::mulTrivial(const BigInt& o) const {
    if (isNaN() || isZero()) return *this;
    if (o.isNaN() || o.isZero()) return o;

    BigInt res;

    uint_long s = _rep.size(), os = o._rep.size(), sos = s + os;
    res.setSize(sos);

    // 外のループを短くしたい
    const BigInt* a;
    const BigInt* b;
    if (s < os) {
        a = this;
        b = &o;
    } else {
        a = &o;
        b = this;
        s = os;
    }

    for (uint_long i = 0; i < s; i++) {
        res += BigInt(*b).operator*=(a->_rep[i]).ushift(i);
    }

    res._sign = _sign != o._sign;

    return res;
}
/*
 * }}}
 */

/*
 * Division & Remainder {{{
 */

// 除算
inline BigInt& BigInt::operator/=(const BigInt& o) { return (*this = div(o)); }
inline BigInt BigInt::div(const BigInt& o, BigInt* r) const {
    BigInt q, rem;
    if (r == NULL) r = &rem;
    divrem(*this, o, &q, r);
    return q;
}
// 剰余
inline BigInt& BigInt::operator%=(const BigInt& o) { return (*this = rem(o)); }
inline BigInt BigInt::rem(const BigInt& o, BigInt* q) const {
    BigInt quot, r;
    if (q == NULL) q = &quot;
    divrem(*this, o, q, &r);
    return r;
}
// 除算剰余
inline bool BigInt::divrem(const BigInt& numer, const BigInt& denom, BigInt* q,
                           BigInt* r) {
    return divremTrivial(numer, denom, q, r);
}
// numerがdenomより絶対値で小さかったらrとしてnumerをそのまま返す
inline bool BigInt::divremTrivial(const BigInt& numer, const BigInt& denom,
                                  BigInt* q, BigInt* r) {
    if (q == NULL || r == NULL) return false;
    if (numer.isNaN() || denom.isNaN() || denom.isZero()) {
        q->NaN();
        r->NaN();
        return false;
    }
    if (numer.compAbs(denom) < 0) {
        *q = ZERO;
        *r = numer;
        return true;
    }

    // u / v
    *r = numer;
    BigInt& u = *r;
    BigInt v = denom;
    auto&& _u = u._rep;
    auto&& _v = v._rep;

    u.normalize();
    v.normalize();

    // 符号を記憶
    bool usign = u._sign;
    bool vsign = v._sign;
    u._sign = v._sign = false;

    // std::cout << "params" << std::endl;
    // std::cout << u.toStrAsVector() << std::endl;
    // std::cout << v.toStrAsVector() << std::endl;

    // u: (m + n)-place, v: n-place
    // n,m >= 1
    const uint_long n = _v.size();
    const uint_long m = _u.size() - n;

    q->setSize(m + 1);

    // v[n-1]が2^(BIT_SIZE-1)以上になるように正規化
    uint_long d = BigInt::numOfLeadingZeros(_v[n - 1]);
    u <<= d;
    v <<= d;
    // u._rep.size()をm + n + 1に
    if (_u.size() == m + n) _u.push_back(0);

    uint_long v1 = _v[n - 1];
    uint_long v0 = n > 1 ? _v[n - 2] : 0;

    v.ushift(m);

    // std::cout << "shifted params" << std::endl;
    // std::cout << u.toStrAsVector() << std::endl;
    // std::cout << v.toStrAsVector() << std::endl;
    // std::cout << m << std::endl;
    // std::cout << n << std::endl;
    // std::cout << d << std::endl;

    for (uint_rep j = m + n; j >= n; --j) {

        if (_u.size() > j) {

            // guess quotient

            uint_long u2 = _u[j];
            uint_long u1 = j > 0 ? _u[j - 1] : 0;
            uint_long u0 = j > 1 ? _u[j - 2] : 0;

            uint_long uhat = u2 << BIT_SIZE | u1;  // 64bit
            uint_long qhat = uhat / v1;            // 32bit or 2^32
            uint_long rhat = uhat - qhat * v1;     // 32bit

            // std::cout << "qhat" << std::endl;
            // std::cout << qhat << std::endl;

            while (!(rhat & UPPER_MASK) &&
                   (qhat & UPPER_MASK ||
                    (qhat * v0 > ((rhat << BIT_SIZE) | u0)))) {
                --qhat;
                rhat += v1;
            }
            // std::cout << qhat << std::endl;

            u -= v * qhat;

            while (u.isNegative()) {
                u += v;
                --qhat;
            }
            // std::cout << qhat << std::endl;

            q->_rep[j - n] = qhat;
        }

        v.lshift();
    }

    u >>= d;

    // 符号を適用
    q->_sign = usign != vsign;
    // 剰余の符号は剰余の符号は被除数準拠
    r->_sign = usign;
    //  u =   v *  q  + r (q >= 0, 0 <= r < v)
    //  u = (-v)*(-q) + r
    // -u =   v *(-q) - r
    // -u = (-v)*  q  - r

    return true;
}

/*
 * }}}
 */

/*
 * Greatest Common Divisor {{{
 */

inline BigInt BigInt::gcd(const BigInt& a, const BigInt& b) {
    return binaryGcd(a, b);
}
inline BigInt BigInt::basicGcd(const BigInt& a, const BigInt& b) {
    if (a.isNaN() || b.isNaN()) return BigInt().NaN();

    BigInt u[2] = {a, b};
    uint_long i = 1;
    for (uint_long j = 2; j--;)
        if (u[j].isNegative()) u[j].flip();

    while (!u[i].isZero()) {
        u[i ^ 1] %= u[i];
        i ^= 1;
    }

    return u[i ^ 1];
}
inline BigInt BigInt::binaryGcd(const BigInt& a, const BigInt& b) {
    if (a.isNaN() || b.isNaN()) return BigInt().NaN();

    if (a.isZero()) return b;
    if (b.isZero()) return a;

    uint_long za = a.numOfTrailingZeros();
    uint_long zb = b.numOfTrailingZeros();
    uint_long k = std::min(za, zb);
    BigInt u[2] = {a >> za, b >> zb};

    while (u[0] != u[1]) {
        uint_long i = u[0] > u[1] ? 0 : 1;
        u[i] -= u[i ^ 1];
        u[i] >>= u[i].numOfTrailingZeros();
    }

    return u[0] << k;
}

/*
 * }}}
 */

/*
 * Montgomery Arithmetic {{{
 */

inline BigInt BigInt::mulMod(const BigInt& o, const BigInt& n) {
    return mulMod(*this, o, n);
}
inline BigInt BigInt::mulMod(const BigInt& a, const BigInt& b,
                             const BigInt& n) {
    if (a.isNaN() || b.isNaN() || n.isNaN() || !n.isPositive())
        return BigInt().NaN();

    if (n.isOdd() && (_M.getModulus() == n || _M.setParams(n)))
        return mulMod(a, b, _M);

    return (a * b) % n;
}
inline BigInt BigInt::mulMod(const BigInt& o, const MontgomerySystem& M) {
    return mulMod(*this, o, M);
}
inline BigInt BigInt::mulMod(const BigInt& a, const BigInt& b,
                             const MontgomerySystem& M) {
    if (a.isNaN() || b.isNaN() || !M.isSet()) return BigInt().NaN();
    return M.toME(M.reduce(a * b));
}
inline BigInt BigInt::powMod(uint_long exp, const BigInt& n) {
    return powMod(*this, exp, n);
}
inline BigInt BigInt::powMod(const BigInt& a, uint_long exp, const BigInt& n) {
    if (a.isNaN() || n.isNaN() || !n.isPositive()) return BigInt().NaN();

    if (n.isOdd() && (_M.getModulus() == n || _M.setParams(n)))
        return powMod(a, exp, _M);

    if (a.isZero()) {
        if (exp == 0) return BigInt().NaN();
        return a;
    }

    BigInt seed = a % n;
    BigInt res((exp & 1) ? seed : BigInt::ONE);

    while (exp >>= 1, exp) {
        seed = seed.sqr() % n;
        if (exp & 1) res = res * seed % n;
    }

    return res;
}
inline BigInt BigInt::powMod(uint_long exp, const MontgomerySystem& M) {
    return powMod(*this, exp, M);
}
inline BigInt BigInt::powMod(const BigInt& a, uint_long exp,
                             const MontgomerySystem& M) {
    if (a.isNaN() || !M.isSet()) return BigInt().NaN();

    if (a.isZero()) {
        if (exp == 0) return BigInt().NaN();
        return a;
    }

    BigInt seed = M.toME(a);
    BigInt res((exp & 1) ? seed : BigInt::ONE);

    while (exp >>= 1, exp) {
        seed = seed.mulMod(seed, M);
        if (exp & 1) res = res.mulMod(seed, M);
    }

    return res;
}

/*
 * }}}
 */

/*
 * Other operations {{{
 */

inline BigInt BigInt::pow(uint_long exp) const {
    if (isNaN()) return *this;

    if (isZero()) {
        if (exp == 0) return BigInt().NaN();
        return *this;
    }

    BigInt seed(*this);
    BigInt res((exp & 1) ? seed : ONE);

    while (exp >>= 1, exp) {
        seed *= seed;
        if (exp & 1) res *= seed;
    }

    return res;
}
inline BigInt BigInt::pow(const BigInt& exp) const {
    if (isNaN()) return *this;

    if (isZero()) {
        if (exp.isZero()) return BigInt().NaN();
        return *this;
    }

    BigInt seed(*this);
    BigInt res(exp.isOdd() ? seed : ONE);
    BigInt e = exp;

    while (e >>= 1, !e.isZero()) {
        seed *= seed;
        if (e.isOdd()) res *= seed;
    }

    return res;
}

/*
 * sqrt
 */

// にゅーとんほう
// std::sqrtは53bitの精度しかない
inline uint_rep BigInt::sqrt(uint_long x) {
    if (x == 0) return 0;
    uint_long s = 1UL << (BigInt::bitLength(x) >> 1);
    do {
        s = (s + x / s) >> 1;
    } while (s > x / s || x / (s + 1) >= (s + 1));
    return s;
}
inline BigInt BigInt::sqrt() const { return sqrtNewton(); }
inline BigInt BigInt::sqrtNewton() const {
    if (isNegative()) return BigInt().NaN();
    if (isNaN() || isZero()) return *this;

    BigInt s = ONE << (bitLength() >> 1);
    while (s = (s + operator/(s)) >> 1,
           operator<(s.sqr()) || operator>=((s + ONE).sqr()))
        ;
    return s;
}
inline BigInt BigInt::sqrtTrivial() const {
    if (isNegative()) return BigInt().NaN();
    if (isNaN() || isZero()) return *this;

    BigInt res;

    // s = sqrt(u)
    BigInt u = *this;
    auto&& _u = u._rep;

    u.normalize();

    // u: m-place
    // n,m >= 1
    uint_long m = _u.size();

    // _u[2n-1] >= floor(b/4)を保証
    uint_long d = BIT_SIZE - BigInt::bitLength(_u[m - 1]);
    if (m & 1) {
        d += BIT_SIZE;
        ++m;
    }
    if (d & 1) d -= 1;
    // u: 2n-place (n >= 1)
    u <<= d;
    uint_long dHalf = d >> 1;

    // init (m = 2n (n >= 1))
    // const uint_long n = m >> 1;
    // uint_long j = n - 1;
    uint_long j = m - 2;

    uint_long _uhead = ((uint_long)_u[j + 1] << BIT_SIZE) | _u[j];  // 64bit
    BigInt s(BigInt::sqrt(_uhead));                                 // 32bit
    BigInt r = BigInt(_uhead) - s.sqr();

    while (j) {
        j -= 2;

        BigInt quot, rem;
        divrem(r.pushLower(_u[j + 1]), s << (dHalf + 1), &quot,
               &rem);  // r broken

        r = rem.pushLower(_u[j]) - quot.sqr();  // rem broken
        while (r.isNegative()) {
            r += (quot << 1) - BigInt::ONE;
            quot -= BigInt::ONE;
        }

        // ありえないけどね！！！
        // if (quot.isNegative() || quot._rep.size() > 1) error

        s.pushLower(quot._rep[0]);
    }

    return s >>= dHalf;
}
inline BigInt BigInt::factorial() const {
    if (isNaN()) return *this;
    if (isNegative()) return ZERO;

    // 中央値
    BigInt center = (*this >> 1) + ONE;
    // 階乗
    BigInt fact = center;
    // カウンタ
    BigInt n = center - (isOdd() ? ONE : TWO);

    center *= center;

    BigInt delta = ONE;
    while (n.isPositive()) {
        center -= delta;
        fact *= center;
        delta += TWO;
        n -= ONE;
    }

    return fact;
}

/*
 * }}}
 */

/*
 * Radix Conversion {{{
 */

inline BigInt BigInt::parse(const std::string& str, uint_long radix) {
    BigInt res;
    parse(res, str, radix);
    return res;
}
inline bool BigInt::parse(BigInt& bint, std::string str, uint_long radix) {
    // ゼロ初期化
    bint.Zero();

    // varidation check
    uint_long size = str.size();
    if (size == 0) return true;
    bool sign = false;
    if (str[0] == '-') {
        sign = true;
        str.erase(0, 1);
        if (--size == 0) {
            bint.NaN();
            return false;
        }
    }

    if (radix < 2 || NUMBER_LIST.size() < radix) {
        bint.NaN();
        return false;
    }
    // to upper case
    std::transform(begin(str), end(str), begin(str), ::toupper);

    if (str.find_first_not_of(NUMBER_LIST.substr(0, radix)) !=
        std::string::npos) {
        bint.NaN();
        return false;
    }

    // determine bases
    uint_long baseSize;
    BigInt base;
    if (radix == 10) {
        baseSize = DEC_SIZE;
        base = BASE10;
    } else {
        baseSize = BIT_SIZE * (std::log(2) / std::log(radix));
        if (baseSize == 0) {
            bint.NaN();
            return false;
        }
        base = BigInt(radix).pow(baseSize);
    }

    // main routine
    uint_long q = size / baseSize, r = size % baseSize;

    if (r) {
        bint += BigInt((int_long)std::stoul(str.substr(0, r), 0, radix));
    }

    for (uint_long i = 0; i < q; ++i) {
        bint *= base;
        bint += BigInt((int_long)std::stoul(
            str.substr(i * baseSize + r, baseSize), 0, radix));
    }

    bint._sign = sign;

    return true;
}

inline std::string BigInt::uintToStr(uint_long num, uint_long radix) {

    if (radix < 2 || NUMBER_LIST.size() < radix) {
        return "*** Invalid radix ***";
    }
    if (num == 0) return "0";

    std::string str;

    while (num) {
        str += NUMBER_LIST[num % radix];
        num /= radix;
    }
    std::reverse(begin(str), end(str));

    return str;
}
inline std::string BigInt::toStr(uint_long radix) const {
    return toStr(*this, radix);
}
// trivial
inline std::string BigInt::toStr(const BigInt& bint, uint_long radix) {
    if (bint.isNaN()) return "NaN";

    if (radix < 2 || NUMBER_LIST.size() < radix) {
        return "*** Invalid radix ***";
    }

    uint_long baseSize;
    BigInt base;
    if (radix == 10) {
        baseSize = DEC_SIZE;
        base = BASE10;
    } else {
        baseSize = BIT_SIZE * (std::log(2) / std::log(radix));
        if (baseSize == 0) return "*** Invalid radix ***";
        base = BigInt(radix).pow(baseSize);
    }

    bool sign = bint._sign;
    std::string str;
    BigInt numer = bint;
    numer._sign = false;
    BigInt quot, rem;

    while (!numer.isZero()) {
        divrem(numer, base, &quot, &rem);

        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(baseSize)
            << uintToStr(rem._rep[0], radix);
        str.insert(0, oss.str());

        numer = quot;
    }

    uint_long head = 0;
    while (head < str.size() && str[head] == '0') ++head;
    if (head == str.size()) return "0";

    // std::cout << bint.toStrAsVector() << std::endl;

    return (sign ? "-" : "") + str.substr(head);
}
inline std::string BigInt::toStrAsVector() const {
    if (isNaN()) return "NaN";

    BigInt rev(*this);
    rev.normalize();
    std::reverse(begin(rev._rep), end(rev._rep));
    bool isFirst = true;
    std::stringstream oss;
    for (auto&& v : rev._rep) {
        if (isFirst) {
            if (rev._sign) oss << "-";
            oss << "[";
            isFirst = false;
        } else {
            oss << ", ";
        }
        oss << v;
    }
    oss << "]";

    return oss.str();
}

// 出力オーバーロード
inline std::ostream& operator<<(std::ostream& os, const BigInt& bint) {
    if (bint.isNaN()) {
        os << "NaN";
        return os;
    }

    // for debug
    // return os << bint.toStrAsVector();
    return os << bint.toStr();
}

/*
 * }}}
 */

/*
 * Random {{{
 */

// 乱数初期化
inline uint_long BigInt::initRandom(uint_long s) {
    if (s == 0) {
        std::random_device rnd;
        s = rnd();
    }
    BigInt::mt32.seed(s);  // 32bit mt
    BigInt::mt64.seed(s);  // 64bit mt
    return 0;
}
inline BigInt BigInt::randomBound(uint_long n, bool crypt) {
    if (n == 0) return ZERO;
    if (crypt) return BigInt().NaN();
    return BigInt(std::uniform_int_distribution<uint_long>(0, n - 1)(mt64));
}
inline BigInt BigInt::randomBound(const BigInt& n, bool crypt) {
    if (n.isNaN() || n.isZero()) return ZERO;
    if (crypt) return BigInt().NaN();

    uint_long i = n._rep.size();
    BigInt res;
    res.setSize(i);

    uint_long b = 0;
    while (n._rep[b] == 0) ++b;

    bool sameFlag = true;
    while (--i > b && sameFlag) {
        res._rep[i] =
            std::uniform_int_distribution<uint_rep>(0, n._rep[i])(mt32);
        if (res._rep[i] < n._rep[i]) sameFlag = false;
    }
    if (sameFlag) {
        res._rep[i] =
            std::uniform_int_distribution<uint_rep>(0, n._rep[i] - 1)(mt32);
    } else {
        while (i--) {
            res._rep[i] = mt32();
        }
    }

    return res;
}
inline BigInt BigInt::randomBits(uint_long n, bool crypt) {
    if (n == 0) return ZERO;
    if (crypt) return BigInt().NaN();

    uint_long q = n / BIT_SIZE, r = n % BIT_SIZE;
    BigInt res;
    res.setSize(q);

    for (auto&& v : res._rep) {
        v = mt32();
    }

    return res.pushUpper(
        std::uniform_int_distribution<uint_rep>(0, (1U << r) - 1)(mt32));
}
inline BigInt BigInt::randomLength(uint_long n, bool crypt) {
    if (n == 0) return ZERO;
    if (crypt) return BigInt().NaN();

    BigInt res = randomBits(n - 1, crypt);
    uint_long q = n / BIT_SIZE, r = n % BIT_SIZE;
    res.setSize(q + 1);
    res._rep[q] |= 1U << (r & BIT_SIZE);

    return res;
}

/*
 * }}}
 */

/*
 *
 * }}}
 *
 */

#endif
