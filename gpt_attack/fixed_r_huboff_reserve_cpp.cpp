#include <gmpxx.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using Poly = std::vector<mpz_class>;

int reserve_threads() {
    const char* raw = std::getenv("FIXED_R_RESERVE_THREADS");
    if (raw == nullptr) return 1;
    int value = std::atoi(raw);
    if (value < 1) return 1;
    return std::min(value, 16);
}

void trim(Poly& poly) {
    while (poly.size() > 1 && poly.back() == 0) {
        poly.pop_back();
    }
}

Poly read_poly() {
    std::size_t size;
    std::cin >> size;
    Poly out(size);
    for (std::size_t i = 0; i < size; ++i) {
        std::cin >> out[i];
    }
    trim(out);
    return out;
}

Poly add_poly(const Poly& lhs, const Poly& rhs) {
    Poly out(std::max(lhs.size(), rhs.size()));
    for (std::size_t i = 0; i < lhs.size(); ++i) out[i] += lhs[i];
    for (std::size_t i = 0; i < rhs.size(); ++i) out[i] += rhs[i];
    trim(out);
    return out;
}

Poly sub_poly(const Poly& lhs, const Poly& rhs) {
    Poly out(std::max(lhs.size(), rhs.size()));
    for (std::size_t i = 0; i < lhs.size(); ++i) out[i] += lhs[i];
    for (std::size_t i = 0; i < rhs.size(); ++i) out[i] -= rhs[i];
    trim(out);
    return out;
}

Poly scale_poly(const Poly& poly, const mpz_class& scalar) {
    if (scalar == 0) return Poly{0};
    Poly out(poly.size());
    for (std::size_t i = 0; i < poly.size(); ++i) out[i] = scalar * poly[i];
    trim(out);
    return out;
}

Poly mul_poly(const Poly& lhs_raw, const Poly& rhs_raw) {
    if (lhs_raw.size() == 1 && lhs_raw[0] == 0) return Poly{0};
    if (rhs_raw.size() == 1 && rhs_raw[0] == 0) return Poly{0};

    const Poly* lhs = &lhs_raw;
    const Poly* rhs = &rhs_raw;
    if (lhs->size() > rhs->size()) std::swap(lhs, rhs);

    const std::size_t out_size = lhs->size() + rhs->size() - 1;
    const std::size_t work = lhs->size() * rhs->size();
    int thread_count = reserve_threads();
    if (thread_count <= 1 || work < 500000 || lhs->size() < 4) {
        Poly out(out_size);
        for (std::size_t i = 0; i < lhs->size(); ++i) {
            if ((*lhs)[i] == 0) continue;
            for (std::size_t j = 0; j < rhs->size(); ++j) {
                if ((*rhs)[j] != 0) {
                    mpz_addmul(
                        out[i + j].get_mpz_t(),
                        (*lhs)[i].get_mpz_t(),
                        (*rhs)[j].get_mpz_t()
                    );
                }
            }
        }
        trim(out);
        return out;
    }

    thread_count = std::min<int>(thread_count, static_cast<int>(lhs->size()));
    std::vector<Poly> partials;
    partials.reserve(thread_count);
    for (int tid = 0; tid < thread_count; ++tid) {
        partials.emplace_back(out_size);
    }

    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (int tid = 0; tid < thread_count; ++tid) {
        std::size_t begin = lhs->size() * tid / thread_count;
        std::size_t end = lhs->size() * (tid + 1) / thread_count;
        threads.emplace_back([&, tid, begin, end]() {
            Poly& local = partials[tid];
            for (std::size_t i = begin; i < end; ++i) {
                if ((*lhs)[i] == 0) continue;
                for (std::size_t j = 0; j < rhs->size(); ++j) {
                    if ((*rhs)[j] != 0) {
                        mpz_addmul(
                            local[i + j].get_mpz_t(),
                            (*lhs)[i].get_mpz_t(),
                            (*rhs)[j].get_mpz_t()
                        );
                    }
                }
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }

    Poly out(out_size);
    for (const Poly& local : partials) {
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] += local[i];
        }
    }
    trim(out);
    return out;
}

bool positive_coefficients(const Poly& poly) {
    for (const auto& coeff : poly) {
        if (coeff <= 0) return false;
    }
    return true;
}

int main() {
    int r = 0;
    long long reserve_denom = 0;
    if (!(std::cin >> r >> reserve_denom)) {
        std::cerr << "bad header\n";
        return 2;
    }

    Poly l_poly = read_poly();
    Poly z_poly = read_poly();
    Poly a_poly = read_poly();
    Poly m_poly = read_poly();

    Poly eval_prev2{1};
    Poly deriv_prev2{0};
    Poly path_eval = eval_prev2;
    Poly log_deriv_num = deriv_prev2;

    if (r > 0) {
        Poly eval_prev1 = add_poly(z_poly, l_poly);
        Poly deriv_prev1 = l_poly;
        path_eval = eval_prev1;
        log_deriv_num = deriv_prev1;

        for (int n = 2; n <= r; ++n) {
            if (n % 2 == 0) {
                path_eval = add_poly(eval_prev1, mul_poly(l_poly, eval_prev2));
                log_deriv_num = add_poly(
                    deriv_prev1,
                    mul_poly(l_poly, add_poly(eval_prev2, deriv_prev2))
                );
            } else {
                path_eval = add_poly(
                    mul_poly(z_poly, eval_prev1),
                    mul_poly(l_poly, eval_prev2)
                );
                log_deriv_num = add_poly(
                    mul_poly(z_poly, deriv_prev1),
                    mul_poly(l_poly, add_poly(eval_prev2, deriv_prev2))
                );
            }
            eval_prev2 = std::move(eval_prev1);
            eval_prev1 = path_eval;
            deriv_prev2 = std::move(deriv_prev1);
            deriv_prev1 = log_deriv_num;
        }
    }

    Poly arm_mean_num = scale_poly(mul_poly(sub_poly(a_poly, Poly{1}), l_poly), 2);
    Poly arm_mean_den = add_poly(z_poly, scale_poly(l_poly, 2));
    mpz_class reserve(static_cast<long>(reserve_denom));
    Poly target_den3 = scale_poly(a_poly, 3 * reserve);
    Poly target_num3 = sub_poly(
        scale_poly(mul_poly(a_poly, sub_poly(Poly{4}, scale_poly(m_poly, 3))), reserve),
        Poly{3}
    );

    Poly term1 = mul_poly(mul_poly(arm_mean_num, path_eval), target_den3);
    Poly term2 = mul_poly(mul_poly(log_deriv_num, arm_mean_den), target_den3);
    Poly arm_eval = mul_poly(arm_mean_den, path_eval);
    Poly term3 = mul_poly(target_num3, arm_eval);
    Poly numerator = add_poly(add_poly(term1, term2), term3);
    Poly denominator = mul_poly(arm_eval, target_den3);

    bool num_ok = positive_coefficients(numerator);
    bool den_ok = positive_coefficients(denominator);
    std::cout << "reserve_ok=" << (num_ok && den_ok ? "true" : "false")
              << " numerator_ok=" << (num_ok ? "true" : "false")
              << " denominator_ok=" << (den_ok ? "true" : "false")
              << " numerator_degree=" << (numerator.size() - 1)
              << " denominator_degree=" << (denominator.size() - 1)
              << "\n";
    return (num_ok && den_ok) ? 0 : 1;
}
