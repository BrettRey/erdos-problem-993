#!/usr/bin/env python3
"""
Testing the key lemma: Product of tail-decreasing sequences.

Claim: If (a_k) and (b_k) are both strictly decreasing for k >= K,
then (c_k) = (a * b)_k is also strictly decreasing for k >= K.
"""

def product_ratio(a, b, k):
    """Compute c_{k+1}/c_k where c = a * b (convolution)."""
    # c_k = sum_{i=0}^k a_i * b_{k-i}
    # c_{k+1} = sum_{i=0}^{k+1} a_i * b_{k+1-i}
    
    # This is expensive to compute exactly, let's approximate
    # In the far tail, the maximum terms dominate
    return "complex"


def simple_test():
    """Test with simple decreasing sequences."""
    # a_k = 1/(k+1) for k >= 0
    a = [1.0 / (k+1) for k in range(20)]
    b = [1.0 / (k+1) for k in range(20)]
    
    # Compute convolution
    c = [0] * 39
    for i in range(len(a)):
        for j in range(len(b)):
            if i + j < len(c):
                c[i+j] += a[i] * b[j]
    
    # Check ratios in tail
    print("a ratios:", [a[i+1]/a[i] for i in range(5, 15)])
    print("b ratios:", [b[i+1]/b[i] for i in range(5, 15)])
    print("c ratios:", [c[i+1]/c[i] for i in range(10, 25)])
    
    # Are they decreasing?
    a_ratios = [a[i+1]/a[i] for i in range(5, 15)]
    b_ratios = [b[i+1]/b[i] for i in range(5, 15)]
    c_ratios = [c[i+1]/c[i] for i in range(10, 25)]
    
    print("\nAre a ratios decreasing?", all(a_ratios[i] >= a_ratios[i+1] for i in range(len(a_ratios)-1)))
    print("Are b ratios decreasing?", all(b_ratios[i] >= b_ratios[i+1] for i in range(len(b_ratios)-1)))
    print("Are c ratios decreasing?", all(c_ratios[i] >= c_ratios[i+1] for i in range(len(c_ratios)-1)))


def test_with_independence_polys():
    """Test with actual independence polynomials."""
    sys.path.insert(0, '.')
    from indpoly import independence_poly
    from targeted import make_broom, make_spider
    
    print("=" * 70)
    print("TESTING PRODUCT RATIOS WITH INDEPENDENCE POLYNOMIALS")
    print("=" * 70)
    
    # Get some polynomials
    for s in [10, 20]:
        n, adj = make_broom(3, s)
        poly = independence_poly(n, adj)
        
        # Find tail region (after mode)
        mode = poly.index(max(poly))
        tail_start = mode + 2
        
        if tail_start < len(poly) - 1:
            tail_ratios = [poly[i+1]/poly[i] for i in range(tail_start, len(poly)-1)]
            print(f"\nBroom(3,{s}): mode={mode}, tail from {tail_start}")
            print(f"  Tail ratios: {[f'{r:.3f}' for r in tail_ratios[:5]]}")
            is_dec = all(tail_ratios[i] >= tail_ratios[i+1] for i in range(len(tail_ratios)-1))
            print(f"  Decreasing: {is_dec}")


if __name__ == '__main__':
    simple_test()
    test_with_independence_polys()
