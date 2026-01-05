import sympy as sp
import sys

# Mengaktifkan pencetakan simbol matematika yang rapi (Unicode)
sp.init_printing(use_unicode=True)

# Unicode symbols
λ = sp.Symbol('λ')

def clean_expr(expr_str):
    """Membersihkan ekspresi matematika untuk tampilan yang lebih rapi"""
    result = str(expr_str)
    result = result.replace('lambda', 'λ')
    # Pangkat: **2 -> ², **3 -> ³
    result = result.replace('**2', '²')
    result = result.replace('**3', '³')
    # Perkalian: * -> ×, tapi bukan di dalam **
    result = result.replace('*', '×')
    return result

def judul():
    print("═"*60)
    print("   KALKULATOR NILAI EIGEN MATRIKS 2x2 DAN 3x3")
    print("   Visualisasi Step-by-Step dengan SymPy")
    print("═"*60)

def input_matriks_sympy(n):
    print(f"\nMasukkan elemen matriks {n}x{n} (pisahkan kolom dengan spasi):")
    elemen = []
    for i in range(n):
        while True:
            try:
                row_input = input(f"  Baris ke-{i+1}: ").strip().split()
                if len(row_input) != n:
                    print(f"  ⚠ Harap masukkan tepat {n} angka.")
                    continue
                row_num = [sp.Rational(x) for x in row_input]
                elemen.append(row_num)
                break
            except ValueError:
                print("  ⚠ Masukkan angka valid.")
    
    return sp.Matrix(elemen)

def format_matrix_str(M, symbol_subs=None):
    """Format matriks sebagai string dengan bracket"""
    rows = M.tolist()
    n = len(rows)
    
    if symbol_subs:
        rows = [[clean_expr(elem) for elem in row] for row in rows]
    else:
        rows = [[str(elem) for elem in row] for row in rows]
    
    # Hitung lebar kolom maksimum
    col_widths = []
    for j in range(len(rows[0])):
        max_w = max(len(rows[i][j]) for i in range(n))
        col_widths.append(max_w)
    
    lines = []
    for i, row in enumerate(rows):
        formatted = "  ".join(row[j].center(col_widths[j]) for j in range(len(row)))
        if i == 0:
            lines.append(f"⎡ {formatted} ⎤")
        elif i == n - 1:
            lines.append(f"⎣ {formatted} ⎦")
        else:
            lines.append(f"⎢ {formatted} ⎥")
    
    return "\n".join(lines)

def tampilkan_langkah(matriks):
    n = matriks.shape[0]
    I = sp.eye(n)

    # ══════════════════════════════════════════════════════════════
    # SOAL
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  SOAL")
    print("═"*60)
    print("\n  Tentukan nilai eigen dari matriks A:")
    print()
    print("  A =")
    print("  " + format_matrix_str(matriks).replace("\n", "\n  "))

    # ══════════════════════════════════════════════════════════════
    # LANGKAH 1: Membentuk A - λI
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  LANGKAH 1: Membentuk A - λI")
    print("═"*60)

    print("\n  A - λI = A - λ·I")
    
    matriks_char = matriks - λ * I
    
    print("\n  A - λI =")
    print("  " + format_matrix_str(matriks_char, symbol_subs=True).replace("\n", "\n  "))

    # ══════════════════════════════════════════════════════════════
    # LANGKAH 2: Persamaan Karakteristik
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  LANGKAH 2: Persamaan Karakteristik det(A - λI) = 0")
    print("═"*60)

    print("\n  PK: det(A - λI) = 0")

    # ══════════════════════════════════════════════════════════════
    # LANGKAH 3: Hitung Determinan
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  LANGKAH 3: Menghitung Determinan")
    print("═"*60)

    polinomial = matriks_char.det()
    polinomial_expanded = sp.expand(polinomial)
    polinomial_str = clean_expr(polinomial_expanded)
    
    if n == 2:
        a, b = matriks_char[0, 0], matriks_char[0, 1]
        c, d = matriks_char[1, 0], matriks_char[1, 1]
        
        a_str = clean_expr(a)
        b_str = clean_expr(b)
        c_str = clean_expr(c)
        d_str = clean_expr(d)
        
        print(f"\n  det = (a)(d) - (b)(c)")
        print(f"\n  det = ({a_str})({d_str}) - ({b_str})({c_str})")
        
        ad = sp.expand(a * d)
        bc = sp.expand(b * c)
        ad_str = clean_expr(ad)
        bc_str = clean_expr(bc)
        
        print(f"\n      = {ad_str} - ({bc_str})")
        print(f"\n      = {polinomial_str}")
        
    elif n == 3:
        print("\n  Ekspansi kofaktor baris pertama:")
        
        a11, a12, a13 = matriks_char[0, 0], matriks_char[0, 1], matriks_char[0, 2]
        
        M11 = matriks_char[1:3, 1:3]
        M12 = matriks_char[1:3, [0, 2]]
        M13 = matriks_char[1:3, 0:2]
        
        det_M11 = sp.expand(M11.det())
        det_M12 = sp.expand(M12.det())
        det_M13 = sp.expand(M13.det())
        
        a11_s = clean_expr(a11)
        a12_s = clean_expr(a12)
        a13_s = clean_expr(a13)
        
        print(f"\n  = ({a11_s})·M₁₁ - ({a12_s})·M₁₂ + ({a13_s})·M₁₃")
        
        det_M11_s = clean_expr(det_M11)
        det_M12_s = clean_expr(det_M12)
        det_M13_s = clean_expr(det_M13)
        
        print(f"\n  = ({a11_s})·({det_M11_s}) - ({a12_s})·({det_M12_s}) + ({a13_s})·({det_M13_s})")
        print(f"\n  = {polinomial_str}")

    print(f"\n  Polinomial Karakteristik:")
    print(f"\n  {polinomial_str} = 0")

    # ══════════════════════════════════════════════════════════════
    # LANGKAH 4: Faktorisasi
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  LANGKAH 4: Faktorisasi Polinomial")
    print("═"*60)

    polinomial_factored = sp.factor(polinomial)
    factored_str = clean_expr(polinomial_factored)
    
    print(f"\n  {polinomial_str} = 0")
    print(f"\n  Difaktorkan:")
    print(f"\n  {factored_str} = 0")

    # ══════════════════════════════════════════════════════════════
    # LANGKAH 5: Nilai Eigen
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  LANGKAH 5: Menentukan Nilai Eigen")
    print("═"*60)

    nilai_eigen = sp.solve(polinomial, λ)
    
    print("\n  Dari faktorisasi, diperoleh:")
    
    for i, val in enumerate(nilai_eigen, 1):
        if val >= 0:
            print(f"\n  λ - {val} = 0  →  λ₍{i}₎ = {val}")
        else:
            print(f"\n  λ + {abs(val)} = 0  →  λ₍{i}₎ = {val}")

    # ══════════════════════════════════════════════════════════════
    # KESIMPULAN
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  KESIMPULAN")
    print("═"*60)
    
    eigen_str = ", ".join([str(v) for v in nilai_eigen])
    print(f"\n  Nilai Eigen = {{ {eigen_str} }}")
    
    print("\n  Dapat ditulis:")
    subscripts = ['₁', '₂', '₃']
    for i, val in enumerate(nilai_eigen):
        sub = subscripts[i] if i < len(subscripts) else f"_{i+1}"
        print(f"    λ{sub} = {val}")

    # ══════════════════════════════════════════════════════════════
    # LANGKAH 6: Vektor Eigen
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  LANGKAH 6: Menghitung Vektor Eigen")
    print("  Untuk setiap λ, selesaikan (A - λI)v = 0")
    print("═"*60)

    subscripts = ['₁', '₂', '₃']
    var_names = ['x', 'y', 'z'] if n == 3 else ['x', 'y']
    
    for idx, val in enumerate(nilai_eigen):
        sub = subscripts[idx] if idx < len(subscripts) else f"_{idx+1}"
        print(f"\n  {'─'*56}")
        print(f"  Untuk λ{sub} = {val}")
        print(f"  {'─'*56}")
        
        # Substitusi nilai eigen
        matriks_substitusi = matriks_char.subs(λ, val)
        
        print(f"\n  [a] Substitusi λ = {val} ke (A - λI):")
        print("  " + format_matrix_str(matriks_substitusi).replace("\n", "\n  "))
        
        # Menampilkan sistem persamaan
        print(f"\n  [b] Sistem Persamaan Linear (A - λI)v = 0:")
        for i in range(n):
            terms = []
            for j in range(n):
                coef = matriks_substitusi[i, j]
                if coef != 0:
                    if coef == 1:
                        terms.append(f"{var_names[j]}")
                    elif coef == -1:
                        terms.append(f"-{var_names[j]}")
                    else:
                        terms.append(f"({coef}){var_names[j]}")
            if terms:
                eq = " + ".join(terms).replace("+ -", "- ")
                print(f"      {eq} = 0")
            else:
                print(f"      0 = 0")
        
        # RREF (Row Reduced Echelon Form)
        rref_matrix, pivots = matriks_substitusi.rref()
        
        print(f"\n  [c] Bentuk Eselon Baris Tereduksi (RREF):")
        print("  " + format_matrix_str(rref_matrix).replace("\n", "\n  "))
        
        # Mencari Null Space
        vektor_basis = matriks_substitusi.nullspace()
        
        print(f"\n  [d] Menentukan Vektor Eigen:")
        
        if len(vektor_basis) > 0:
            # Tentukan variabel bebas dan terikat
            free_vars = [j for j in range(n) if j not in pivots]
            
            if free_vars:
                free_var_names = [var_names[j] for j in free_vars]
                print(f"      Variabel bebas: {', '.join(free_var_names)}")
                print(f"      Misalkan {free_var_names[0]} = t (parameter)")
            
            for vec in vektor_basis:
                print(f"\n      Solusi umum:")
                # Tampilkan vektor dalam bentuk parameter
                vec_elems = [str(vec[i]) for i in range(vec.rows)]
                
                # Tampilkan sebagai vektor kolom dengan keterangan
                print(f"      v{sub} = t ×", end="")
                
                # Normalisasi vektor (ambil elemen terakhir sebagai referensi)
                for i, elem in enumerate(vec_elems):
                    if i == 0:
                        print(f" ⎡ {elem} ⎤")
                    elif i == vec.rows - 1:
                        print(f"             ⎣ {elem} ⎦")
                    else:
                        print(f"             ⎢ {elem} ⎥")
            
            print(f"\n  [e] Basis Ruang Eigen untuk λ{sub} = {val}:")
            for vec in vektor_basis:
                print(f"\n      v{sub} =")
                print("      " + format_matrix_str(vec).replace("\n", "\n      "))
        else:
            print("      Tidak ditemukan vektor non-trivial.")
            print("      (Matriks mungkin defektif)")

    # ══════════════════════════════════════════════════════════════
    # KESIMPULAN AKHIR
    # ══════════════════════════════════════════════════════════════
    print("\n" + "═"*60)
    print("  RINGKASAN HASIL")
    print("═"*60)
    
    print("\n  Nilai Eigen dan Vektor Eigen:")
    for idx, val in enumerate(nilai_eigen):
        sub = subscripts[idx] if idx < len(subscripts) else f"_{idx+1}"
        matriks_substitusi = matriks_char.subs(λ, val)
        vektor_basis = matriks_substitusi.nullspace()
        
        print(f"\n  λ{sub} = {val}")
        if vektor_basis:
            for vec in vektor_basis:
                vec_str = ", ".join([str(vec[i]) for i in range(vec.rows)])
                print(f"  v{sub} = ({vec_str})ᵀ")

    print("\n" + "═"*60)
    print("  SELESAI")
    print("═"*60)

def main():
    while True:
        judul()
        print("\nPilih Ordo Matriks:")
        print("  1. Matriks 2×2")
        print("  2. Matriks 3×3")
        print("  3. Keluar")
        
        pilihan = input("\nPilihan (1/2/3): ")

        if pilihan == '1':
            M = input_matriks_sympy(2)
            tampilkan_langkah(M)
        elif pilihan == '2':
            M = input_matriks_sympy(3)
            tampilkan_langkah(M)
        elif pilihan == '3':
            sys.exit()
        else:
            print("  Pilihan tidak valid.")
        
        input("\nTekan Enter untuk lanjut...")
        print("\n" * 2)

if __name__ == "__main__":
    main()