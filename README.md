# Linear Programming Solver - Giải bài toán Quy hoạch tuyến tính tổng quát

Chương trình C++ giải bài toán **Quy hoạch tuyến tính (QHTT)** tổng quát qua giao diện dòng lệnh, hỗ trợ ba thuật toán kinh điển và in chi tiết **từng bước lặp của từ vựng** (bảng đơn hình) để người học dễ dàng đối chiếu với lời giải tay.

## ✨ Tính năng

- **Ba phương pháp giải:**

  | # | Phương pháp | Điều kiện áp dụng | Đặc điểm |
  |---|---|---|---|
  | 1 | Đơn hình (Simplex) | Tất cả $b_i \ge 0$ | Nhanh, kinh điển |
  | 2 | Quy tắc Bland | Tất cả $b_i \ge 0$ | Tránh lặp vô hạn (cycling) khi suy biến |
  | 3 | Hai pha (Two-phase) | Mọi trường hợp | Tổng quát nhất |

  Có thể chọn chạy **tất cả** phương pháp áp dụng được (tùy chọn 4); chương trình tự kiểm tra tính hợp lệ của phương pháp đã chọn.

- **Xử lý bài toán tổng quát:** Min/Max, ràng buộc ≥ / ≤ / =, biến không âm, biến không dương, biến tự do.
- **Phát hiện đầy đủ các trường hợp:** nghiệm tối ưu duy nhất, vô số nghiệm (in biểu thức tham số theo biến tự do kèm điều kiện), vô nghiệm, không giới nội, lặp vô hạn.
- **Minh họa tường minh:** in từng từ vựng qua mỗi vòng lặp, biến vào, biến ra, tỉ lệ tại các ràng buộc như giải tay.
- **Kiểm tra nhập liệu:** yêu cầu nhập lại khi sai định dạng.

## 🛠️ Biên dịch và chạy

```bash
g++ -std=c++17 -o lp_solver linear_programming_solver.cpp
./lp_solver
```

Trên Windows (MinGW): `g++ -std=c++17 -o lp_solver.exe linear_programming_solver.cpp`

## 📋 Luồng nhập liệu

1. **Loại tối ưu:** `1` = Min, `2` = Max
2. **Số biến** và **số ràng buộc**
3. **Hệ số hàm mục tiêu** (lần lượt từng biến)
4. **Mỗi ràng buộc:** hệ số từng biến → vế phải $b_i$ → dấu (`1`: ≥, `2`: ≤, `3`: =)
5. **Điều kiện mỗi biến:** `1`: $x_i \ge 0$, `2`: $x_i \le 0$, `3`: tự do
6. **Phương pháp giải:** `1` Đơn hình, `2` Bland, `3` Hai pha, `4` Tất cả

## 📄 Tài liệu

- **`Bao_cao_Quy_hoach_tuyen_tinh.pdf`** - báo cáo đầy đủ: hướng dẫn sử dụng chi tiết kèm ví dụ minh họa, bảng tra cứu ký hiệu, xử lý lỗi thường gặp, ưu/khuyết điểm và so sánh với các phần mềm QHTT khác (LINGO, CPLEX, Gurobi, GLPK, MATLAB, OR-Tools).
- Phương pháp trình bày hai pha dạng bảng tham khảo tại [bảng tính minh họa](https://docs.google.com/spreadsheets/d/1IVz3A70W21_D2cxX5cuoR7_7Pmy0Nt3D_tIZdt1fDwY/edit?usp=sharing).

## ⚠️ Hạn chế hiện tại

- Nhập liệu thủ công bằng văn bản, dễ nhầm với bài toán nhiều biến; nhập sai phải chạy lại từ đầu.
- Dùng con trỏ thô và cấp phát động (`new`/`delete[]`) ở phần nhập liệu.
- Cài đặt phục vụ học tập (bài toán cỡ sách giáo khoa), chưa tối ưu hiệu năng.
- Trường hợp vô số nghiệm với 2+ biến tự do: điều kiện biến tự do chỉ in ở dạng đơn giản (≥ 0).
