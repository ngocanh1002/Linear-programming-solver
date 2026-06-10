#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <float.h>
#include <cmath>
#include <algorithm>
#include <climits>

/* Định nghĩa các biến
1. Biến Đầu Vào & Mô Hình Tổng Quát
    - a/check: Loại tối ưu: 1: Min, 2: Max
    - b: số lượng biến chính ban đầu (x1,x2,...)
    - c: số lượng ràng buộc ban đầu
    - a_heso[i]: hệ số hàm mục tiêu của biến chính xi (cj ban đầu)
    - rangbuoc_heso[i][j]: mảng 2 chiều chứa hệ số các biến trong từng ràng buộc
    - sothucbi[i]: Mảng chứa các giá trị vế phải của ràng buộc (vector số thực bi ban đầu).
    - rangbuoc_dau[i]: dấu của ràng buộc thứ i (1: ≥, 2: ≤, 3: =)
    - dieukien_bien[i]: dấu của từng biến chính (xi) riêng lẻ (1: ≥, 2: ≤, 3: tự do)
2. Biến trong hàm MoHinhChuanTac
    - loai_toi_uu: lưu lại loại bài toán (giống check)
    - so_bien_goc , so_rang_buoc_goc: lưu lại số lượng biến và ràng buộc ban đầu để tham chiếu
    - so_bien_moi, so_rang_buoc_moi: lưu lại số lượng biến và ràng buộc mới sau chuẩn hóa (vd: biến tự do xj thành xj+ và xj-)
    - ten_bien: vector string lưu tên các biến mới sau chuẩn hóa ("x1", "x2", "y1", "y2", "x1+", "x1-", "x2+", "x2-", ...)
    - heso_ham_muc_tieu: vectorlưu hệ số của hàm mục tiêu đã chuẩn hóa (đổi dấu nếu là bài toán Max)
    - heso_rang_buoc: mảng 2 chiều chứa hệ số các ràng buộc đã chuẩn hóa
    - b_sothuc: vector số thực b đã chuẩn hóa (đối dấu cả vế ràng buộc nếu "≥")
    - lich_su_kiemtra, da_kiemtra: Các biến mutable dùng để lưu kết quả kiểm tra khả năng áp dụng các phương pháp, giúp không phải tính lại khi chọn "Tất cả phương pháp".
3. Biến trong các Hàm Giải Thuật (Đơn hình & Hai pha)
    - tonghop, tonghop_pha2: Bảng đơn hình (ma trận 2 chiều) chứa toàn bộ hệ số của bài toán, bao gồm cả các biến phụ và cột vế phải.
    - indices_bien_co_so: lưu chỉ số (index) các biến đang là biến cơ sở, tức là nằm trên trục dọc trong bảng đơn hình / biến phụ
    - cb: Vector lưu hệ số trong hàm mục tiêu của các biến đang nằm trong cơ sở. Dùng để tính hàng zj
    - heso_hammuctieu, original_hammuctieu: Vector lưu hệ số của hàm mục tiêu z (bao gồm biến chính và biến phụ)
    - z: Vector tạm để lưu các giá trị của hàng zj được tính ở mỗi vòng lặp
    - hammuctieu (biến cục bộ): Vector tạm để lưu hàng mục tiêu cj − zj . Dựa vào đây để xác định tính tối ưu và tìm biến vào
    - pivot_col, pivot_row: chỉ số cột(biến vào) và hàng(biến ra) của phần tử xoay (pivot key) trong bảng đơn hình
    - bien_tu_do: Vector lưu chỉ số của các biến phi cơ sở có hệ số bằng 0 trong hàng cj−zj , dùng để xác định bài toán có vô số nghiệm.
4. Biến trong hàm xu_ly_luachon
    - phphap: lưu lại lựa chọn phương pháp giải (1: đơn hình, 2: bland, 3: hai pha)
    - phphap_chon: lưu lại lựa chọn phương pháp giải (1: đơn hình, 2: bland, 3: hai pha)
5. Biến trong hàm kiemtra_phuongphap
    - ketqua: Vector lưu kết quả kiểm tra phương pháp (0: không áp dụng, 1: áp dụng)
    - tontai_bi_am: kiểm tra xem có tồn tại số thực bi âm không
6. Biến trong hàm kiem_tra_trung_lap
    - current_table: Bảng đơn hình hiện tại
    - current_basis: Vector lưu chỉ số các biến cơ sở hiện tại
    - initial_table: Bảng đơn hình ban đầu
7. Biến trong hàm giai_phuongphap_tieuchuan
    - mo_hinh: mô hình chuẩn tắc
    - use_bland: kiểm tra xem phương pháp Bland đã được chọn hay chưa
    - use_2pha: nếu như tất cả số thực bi đều >0 mà use 2 pha thì chuyển sang pha 2 (base là pp đơn hình)
*/

using namespace std;

const double epsilon = 1e-9; // sai số cho phép

// Cấu trúc lưu trữ dữ liệu mô hình chuẩn tắc
struct MoHinhChuanTac{
    int loai_toi_uu; // 1: Min, 2:Max
    int so_bien_goc; // tổng số biến gốc (xi)
    int so_bien_moi;  // tổng số biến mới (xi,yi)
    int so_rang_buoc_goc; // tổng số ràng buộc gốc
    int so_rang_buoc_moi; // tổng số ràng buộc mới
    vector<double> heso_ham_muc_tieu; // hệ số của hàm mục tiêu (z) sau chuẩn hóa
    vector<vector<double>> heso_rang_buoc; // hệ số của các ràng buộc sau chuẩn hóa
    vector<double> b_sothuc; // vector bên vế phải của ràng buộc
    vector<string> ten_bien; // tên biến sau chuẩn hóa (x,y)
    // mutable: cho phép thay đổi giá trị của biến trong hàm const (hàm const: hàm k thể thay đổi dữ liệu của đối tượng)
    mutable vector <bool> lich_su_kiemtra; // lưu lại kết quả kiểm tra phương pháp để không phải kiểm tra lại cho chọn option 4. All phương pháp
    mutable bool da_kiemtra = false; // đánh dấu đã kiểm tra phương pháp chưa
    
    // constructor
    MoHinhChuanTac(): loai_toi_uu(1), so_bien_goc(0), so_bien_moi(0), so_rang_buoc_goc(0), so_rang_buoc_moi(0) {};
    
    // Xuất ra các dữ liệu của mô hình đã chuẩn hóa
    void in_du_lieu(){
        cout<<"\n\nDữ liệu của mô hình chuẩn tắc khi đã chuẩn hóa:\n";
        cout<<"Tổng số biến mới: "<<so_bien_moi<<endl;
        cout<<"Tổng số ràng buộc mới: "<<so_rang_buoc_moi<<endl;
        cout<<"Tên các biến sau chuẩn hóa (x^T): [ ";
            for(int i=0;i<ten_bien.size();i++){
                cout<<ten_bien[i]<<" ";
            }
            cout <<" ]"<<endl;
        
        cout<<"Hệ số hàm mục tiêu (c^T): [ ";
            for(int i=0; i<heso_ham_muc_tieu.size(); i++){
                cout<<heso_ham_muc_tieu[i]<<" ";
            }
            cout<<" ]"<<endl;
            
        cout<<"Ma trận hệ số A: \n[ ";
        for(int i=0; i<heso_rang_buoc.size(); i++){
            for(int j=0; j<heso_rang_buoc[i].size(); j++){
                cout<<heso_rang_buoc[i][j] <<" ";
            }
            if(i<heso_rang_buoc.size()-1) {
                cout<<endl;
                cout<<"  ";
            }
        }
            cout<<" ]"<<endl;
        
        cout<<"Vector b: [ ";
            for(int i=0; i<b_sothuc.size(); i++){
                cout<<b_sothuc[i]<<" ";
            }
            cout<<" ]" << endl;
            cout<<"\n";
    }
};

// Hàm hiển thị mô hình tổng quát
void mohinh_tongquat(int check, int b, int c, double* a_heso,double** rangbuoc_heso, int* rangbuoc_dau, double* sothucbi,int* dieukien_bien){
    // Hiển thị lại hệ số hàm mục tiêu và các hệ số ràng buộc
    cout << "\n============================\n";
    cout <<  "║     Mô hình tổng quát    ║";
    cout << "\n============================\n";
    if(check == 1) cout << "Min  ";
    else cout << "Max  ";

    // In hàm mục tiêu theo định dạng: 4x1 + 4x2 + ...
    for (int i = 0; i < b; i++) {
        cout << a_heso[i] << "x" << i+1;
        if (i < b - 1) {
            cout << " + ";  // Thêm dấu "+" giữa các hệ số, trừ cho phần cuối
        }
    }
    cout << "\n";

    // In các ràng buộc
    for (int j = 0; j < c; j++) {
        for (int i = 0; i < b; i++) {
            if(i==0){
                cout<<"     ";
            }
            cout << rangbuoc_heso[j][i] << "x" << i+1;
            if (i < b - 1) {
                cout << " + ";  // Thêm dấu "+" giữa các hệ số, trừ cho phần cuối
            }
        }
        if(rangbuoc_dau[j]==1)cout<<" ≥ ";
        else if(rangbuoc_dau[j]==2)cout <<" ≤ ";
        else cout<<" = ";
        cout << sothucbi[j];
        cout << endl;
    }
    
    // điều kiện ràng buộc mặc định
    for(int i=0;i<b;i++){
        cout<<"     ";
        cout<<"x"<<i+1;
        if(dieukien_bien[i]==1) cout<<" ≥ 0";
        else if(dieukien_bien[i]==2) cout<<" ≤ 0";
        else cout<< " tự do.";
        cout<<"\n";
    }
    
}

// Hàm chuyển đổi mô hình tổng quát sang mô hình chuẩn tắc, trả về dữ liệu (để use cho các pp giải)
MoHinhChuanTac chuyensang_mohinhchuantac(int check, int b, int c, double* a_heso,double** rangbuoc_heso, int* rangbuoc_dau, double* sothucbi,int* dieukien_bien){
    MoHinhChuanTac ket_qua;
    
    // thông tin cơ bản
    ket_qua.loai_toi_uu = check;
    ket_qua.so_bien_goc = b;
    ket_qua.so_rang_buoc_goc = c;
    
    // so_bien_moi
    int bienmoi=0;
    for(int i=0;i<b;i++){
        if(dieukien_bien[i]==3) bienmoi+=2; // tự do có xj+ và xj-
        else bienmoi+=1;
    }
    ket_qua.so_bien_moi = bienmoi;
    
    //so_rang_buoc_moi
    int rangbuocmoi=c;
    for(int i=0; i<c; i++){
        if(rangbuoc_dau[i]==3) rangbuocmoi++; // nếu là '=' thì chia ra thành <= và >=
    }
    ket_qua.so_rang_buoc_moi = rangbuocmoi;
    
    // tạo các tên biến mới
    ket_qua.ten_bien.resize(bienmoi); // chỉnh sửa kích thước của tên biến theo bienmoi
    int idx_bien=0;
    
    for(int i=0;i<b;i++){
        if(dieukien_bien[i]==3){
            ket_qua.ten_bien[idx_bien] = "x"+ to_string(i+1) + "+"; // phải xài to_string để chuyển từ dạng number thành string vì định dạng của ten_bien là string. 
            idx_bien++;
            ket_qua.ten_bien[idx_bien] = "x" +to_string(i+1) + "-";
            idx_bien++;
        }
        // <= : chuyển thành y và dấu >=
        else if(dieukien_bien[i]==2){
            ket_qua.ten_bien[idx_bien] = "y" + to_string(i+1);
            idx_bien++;
        }
        // >= giữ nguyên
        else {
            ket_qua.ten_bien[idx_bien] = "x" + to_string(i+1);
            idx_bien++;
    }
    }

    // hệ số hàm mục tiêu (c)
    ket_qua.heso_ham_muc_tieu.resize(bienmoi);
    idx_bien=0;
    for(int i=0; i<b; i++){
        double heso = (check==2) ? -a_heso[i] : a_heso[i];
        // TH có xi+ và xi-
        if(dieukien_bien[i]==3){
            ket_qua.heso_ham_muc_tieu[idx_bien] = heso;
            idx_bien++;
            ket_qua.heso_ham_muc_tieu[idx_bien] = -heso; // - (xi-)
            idx_bien++;
        }
        else if(dieukien_bien[i]==2){
            ket_qua.heso_ham_muc_tieu[idx_bien] = -heso;
            idx_bien++;
        }
        else {
            ket_qua.heso_ham_muc_tieu[idx_bien] = heso;
            idx_bien++;
        }
    }
    
    // ma trận các hệ số ràng buộc
    ket_qua.heso_rang_buoc.resize(rangbuocmoi, vector<double>(bienmoi,0.0)); // khởi tạo ma trận với các phần tử mặc định là 0.0
    
    // ma trận bi (vế phải)
    ket_qua.b_sothuc.resize(rangbuocmoi);
    
    int idx_rang_buoc =0;
    for(int i=0; i<c; i++){
        // vế phải _ bi
        if(rangbuoc_dau[i]==1) {
            ket_qua.b_sothuc[idx_rang_buoc] = -sothucbi[i];
        } else{
            ket_qua.b_sothuc[idx_rang_buoc] = sothucbi[i];
        }
    
        // vế trái
        idx_bien = 0;
        for (int j=0; j<b; j++){
            double heso;
            
            if(rangbuoc_dau[i]==1){ // >= đảo dấu <=
                heso= -rangbuoc_heso[i][j];
            } else {
                heso = rangbuoc_heso[i][j];
            }
            
            if(dieukien_bien[j]==3){
                // biến tự do : xj = xj+ - xj-
                ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien] = heso;      // Hệ số cho xj+
                ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien + 1] = -heso; // Hệ số cho xj-
                idx_bien += 2;
            } else if (dieukien_bien[j] == 2) {
                // Biến <= 0: xj = -yj
                ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien] = -heso; // Hệ số cho yj
                idx_bien += 1;
            } else {
                // Biến >= 0: giữ nguyên
                ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien] = heso;
                idx_bien += 1;
            }
        }
        
        idx_rang_buoc++;
        
        // Nếu là ràng buộc =, tạo thêm một ràng buộc với dấu ngược lại
        if (rangbuoc_dau[i] == 3) {
            ket_qua.b_sothuc[idx_rang_buoc] = -sothucbi[i];
            
            // Xử lý vế trái (hệ số các biến) với dấu ngược lại
            idx_bien = 0;
            for (int m = 0; m < b; m++) {
                double heso = -rangbuoc_heso[i][m]; // Đảo dấu
                
                if (dieukien_bien[m] == 3) {
                    // Biến tự do: xj = xj+ - xj-
                    ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien] = heso;      // Hệ số cho xj+
                    ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien + 1] = -heso; // Hệ số cho xj-
                    idx_bien += 2;
                } else if (dieukien_bien[m] == 2) {
                    // Biến <= 0: xj = -yj
                    ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien] = -heso; // Hệ số cho yj
                    idx_bien += 1;
                } else {
                    // Biến >= 0: giữ nguyên
                    ket_qua.heso_rang_buoc[idx_rang_buoc][idx_bien] = heso;
                    idx_bien += 1;
                }
            }
            
            idx_rang_buoc++;
        }
    }
    
    return ket_qua;
}

// Hàm hiển thị mô hình chuẩn tắc
void hienthi_MoHinhChuanTac(const MoHinhChuanTac& mo_hinh){
    // const: mo_hinh k thể bị thay đổi trong hàm -> tức chỉ use dữ liệu từ mo_hinh mà k sửa đổi nó
    cout << "\n===========================\n";
    cout <<  "║     Mô hình chuẩn tắc   ║";
    cout << "\n===========================\n";
    
    // Hiển thị hàm mục tiêu
    if(mo_hinh.loai_toi_uu==1) cout<<"Min ";
    else cout<<"- Min ";  // Max -> - Min
    
    for(int i=0; i<mo_hinh.so_bien_moi; i++){
        double heso = mo_hinh.heso_ham_muc_tieu[i];
        if (heso!=0){
            if(i==0){ // phần tử in đầu tiên (x1)
                if(heso<0) cout<<"-"<< -heso <<mo_hinh.ten_bien[i];
                else cout << heso<< mo_hinh.ten_bien[i];
            } else {
            if(heso<0) cout << " - "<< -heso <<mo_hinh.ten_bien[i];
            else cout <<" + "<< heso<< mo_hinh.ten_bien[i];
            }
        }
    }
    cout<<endl;
    
    // Hiển thị các ràng buộc
    for(int j=0; j<mo_hinh.so_rang_buoc_moi; j++){
        if(mo_hinh.loai_toi_uu==2) cout<<" ";
        cout << "    ";
        
        for(int i=0; i<mo_hinh.so_bien_moi; i++){
            double heso = mo_hinh.heso_rang_buoc[j][i];
            if(i==0){
                if(heso<0) cout<<"-"<< -heso << mo_hinh.ten_bien[i];
                else cout<<heso << mo_hinh.ten_bien[i];
            } else {
                if(heso<0) cout<<" - "<< -heso << mo_hinh.ten_bien[i];
                else cout<<" + " << heso << mo_hinh.ten_bien[i];
            }
        }
        cout << " ≤ " << mo_hinh.b_sothuc[j] << endl; 
    }
    
    // Hiển thị các ràng buộc không âm cho các biến xi riêng lẻ
    for (int i=0; i<mo_hinh.so_bien_moi; i++){
        if(mo_hinh.loai_toi_uu==2) cout<<" ";
        cout << "    " << mo_hinh.ten_bien[i] << " ≥ 0" << endl;
    }
}

// Hàm tổng hợp kiểm tra khả năng áp dụng các phương pháp
vector <bool> kiemtra_phuongphap(const MoHinhChuanTac& mo_hinh){
    // Nếu đã kiểm tra rồi thì trả về kết quả đã lưu
    if (mo_hinh.da_kiemtra) return mo_hinh.lich_su_kiemtra;

    vector <bool> ketqua(3,false); // 3 phương pháp [0]: đơn hình, [1]: bland, [2]: hai pha.
    // kiểm tra xem có tồn tại số thực bi âm không
    bool tontai_bi_am = any_of(
        mo_hinh.b_sothuc.begin(), // duyệt từ phần tử đầu tiên 
        mo_hinh.b_sothuc.end(),   // đến phần tử cuối cùng
        [] (double b) { // điều kiện kiểm tra: có phần tử nào nhỏ hơn 0 không?
            return b<0;
        }
    );

    ketqua[0] = !tontai_bi_am; // nếu không tồn tại bi âm thì phương pháp đơn hình áp dụng được
    ketqua[1] = !tontai_bi_am; // bland tương tự đơn hình
    ketqua[2] = true; // 2 pha áp dụng được all TH (bi>=0,<=0,=0)

    // Lưu vào lịch sử (để không phải kiểm tra lại cho chọn option 4. All phương pháp)
    mo_hinh.lich_su_kiemtra = ketqua;
    mo_hinh.da_kiemtra = true;

    return ketqua;
}

// Hàm kiểm tra xem có trùng lặp với trạng thái ban đầu ko (từ vựng gốc)
bool kiem_tra_trung_lap(const vector<vector<double>>& current_table,
                       const vector<int>& current_basis,
                       const vector<vector<double>>& initial_table,
                       const vector<int>& initial_basis) { 
    // Kiểm tra số chiều
    if (current_table.size() != initial_table.size() || 
        current_basis.size() != initial_basis.size()) {
        return false;
    }

    // Kiểm tra ma trận hệ số
    for (size_t i = 0; i < current_table.size(); ++i) {
        for (size_t j = 0; j < current_table[i].size(); ++j) {
            if (fabs(current_table[i][j] - initial_table[i][j]) > epsilon) { // vì sự làm tròn của double nên có thể xảy ra sai số
                return false;
            }
        }
    }

    // Kiểm tra biến cơ sở
    for (size_t i = 0; i < current_basis.size(); ++i) {
        if (current_basis[i] != initial_basis[i]) {
            return false;
        }
    }

    return true;
}   

// Hàm giải bài toán QHTT bằng phương pháp đơn hình (biểu diễn dạng bảng)
void giai_phuongphap_tieuchuan(MoHinhChuanTac& mo_hinh, bool use_bland = false, bool use_2pha = false) {
    // =====================================================================
    // [BẢN CÔNG KHAI] Phần cài đặt thuật toán đã được lược bỏ
    // để tránh sao chép. Bản đầy đủ cài đặt tại đây toàn bộ logic:
    //   - Tính z_j và hàng ước lượng (c_j - z_j)
    //   - Chọn biến vào / biến ra (quy tắc Dantzig / Bland)
    //   - Phép xoay (pivot) và cập nhật từ vựng qua từng vòng lặp
    //   - Phát hiện tối ưu / vô số nghiệm / vô nghiệm / không giới nội
    // Liên hệ tác giả để được cung cấp mã nguồn đầy đủ.
    // =====================================================================
    cout << "\n[Bản công khai] Phần giải thuật đã được lược bỏ. "
         << "Vui lòng liên hệ tác giả để xem bản đầy đủ.\n";
}

void giai_phuongphap_2pha(MoHinhChuanTac& mo_hinh) {
    // =====================================================================
    // [BẢN CÔNG KHAI] Phần cài đặt thuật toán đã được lược bỏ
    // để tránh sao chép. Bản đầy đủ cài đặt tại đây toàn bộ logic:
    //   - Tính z_j và hàng ước lượng (c_j - z_j)
    //   - Chọn biến vào / biến ra (quy tắc Dantzig / Bland)
    //   - Phép xoay (pivot) và cập nhật từ vựng qua từng vòng lặp
    //   - Phát hiện tối ưu / vô số nghiệm / vô nghiệm / không giới nội
    // Liên hệ tác giả để được cung cấp mã nguồn đầy đủ.
    // =====================================================================
    cout << "\n[Bản công khai] Phần giải thuật đã được lược bỏ. "
         << "Vui lòng liên hệ tác giả để xem bản đầy đủ.\n";
}

// Hàm xử lý các phương pháp chọn -> để chuyển hướng đến các phương pháp khác nhau
void xu_ly_luachon(MoHinhChuanTac& mo_hinh, int phphap){
    auto kiem_tra = kiemtra_phuongphap(mo_hinh);

    switch (phphap){
        case 1: // đơn hình
            if(kiem_tra[0]) giai_phuongphap_tieuchuan(mo_hinh);
            else cout << "Không thể áp dụng phương pháp đơn hình vì tồn tại số thực b_i < 0\n";
            break;
        case 2: // Bland
            if(kiem_tra[1]) giai_phuongphap_tieuchuan(mo_hinh,true);
            else cout << "Không thể áp dụng phương pháp Bland vì tồn tại số thực b_i < 0\n";
            break;
        case 3: // 2 pha
            giai_phuongphap_2pha(mo_hinh); // vì pp này luôn áp dụng được.
            break;
        case 4: // all phương pháp
            if(kiem_tra[0]) giai_phuongphap_tieuchuan(mo_hinh);
            else cout << "Không thể áp dụng phương pháp đơn hình vì tồn tại số thực b_i < 0\n";
            if(kiem_tra[1]) giai_phuongphap_tieuchuan(mo_hinh,true);
            else cout << "Không thể áp dụng phương pháp Bland vì tồn tại số thực b_i < 0\n";
            
            giai_phuongphap_2pha(mo_hinh);
            break;

        default:
            cout<< "Lựa chọn không hợp lệ~"<<endl;
    }
}


int main() {
    string input;
    double a;
    int b,c;
    double heso;
    
    // Nhập LOẠI TỐI ƯU
    cout << "Nhập loại tối ưu (1 cho min, 2 cho max): " ;
    while(true){
        cin >> input; // nhập vào số or chuỗi
        try{
            a = stod(input);
            if(a==1||a==2) break;
            cout << "Lỗi! Vui lòng nhập số hợp lệ (1 hoặc 2): ";
        }
        catch (invalid_argument& e){
            // Nếu không thể chuyển đổi chuỗi thành số (chắc chắn là chữ)
            cout << "Lỗi! Vui lòng nhập số hợp lệ (1 hoặc 2): ";
        }
    }
    cout <<"\n";
    
    // Nhập số biến
    cout << "Nhập số biến (xi): ";
    while(true){
        cin >> input;
        stringstream ss(input);
        if(ss >> b && ss.eof() && b>0) break;
        cout <<"Lỗi! Vui lòng nhập số nguyên hợp lệ và lớn hơn 0: ";
    }
    cout <<"\n";
    
    // NHẬP SỐ RÀNG BUỘC 
    cout << "Nhập số ràng buộc: ";
    while(true){
        cin >> input;
        stringstream ss(input);
        if(ss >> c && ss.eof()) break;
        cout <<"Lỗi! Vui lòng nhập số nguyên hợp lệ: ";
    }
    cout <<"\n";
    
     // Cấp phát mảng động cho hệ số hàm mục tiêu và các hệ số ràng buộc
    double* a_heso = new double[b];
    double *sothucbi = new double[c];
    double** rangbuoc_heso = new double*[c];
    for (int i=0;i<c;i++){
        rangbuoc_heso[i] = new double[b];
    }
    int* rangbuoc_dau = new int[c];
    int* dieukien_bien = new int[b];
    
    
    // NHẬP HỆ SỐ HÀM MỤC TIÊU
    cout<<"Hàm mục tiêu: \n";
    for(int i=0;i<b;i++){
        cout<<"Hệ số cho x"<<i+1<<": ";
        cin >> a_heso[i];
    }
    cout <<"\n";
    
    // NHẬP HỆ SỐ RÀNG BUỘC
    cout<<"Ràng buộc/ Bất đẳng thức: \n";
    cout<<"---------------------------";
    for(int j=0;j<c;j++){
        cout<<"\nRàng buộc thứ "<<j+1<<": \n";
        for(int i=0;i<b;i++){
            cout<<"Hệ số x"<<i+1<<": ";
            cin >> rangbuoc_heso[j][i];
        }
        cout<<"\nHệ số thực (bi) của ràng buộc: ";
        cin >> sothucbi[j];
    cout <<"\n";
    
    // ≤/=/≥ DẤU RÀNG BUỘC
        cout<<"Dấu của ràng buộc "<<j+1<<"(1: ≥, 2: ≤, 3: =): ";
        while(true){
            cin >> rangbuoc_dau[j];
            if(rangbuoc_dau[j]>=1&&rangbuoc_dau[j]<=3) break;
            cout<<"Lỗi! Vui lòng nhập dấu hợp lệ (1, 2, hoặc 3): ";
        }
    }
    cout<<endl;
    
    // ĐIỀU KIỆN BIẾN RIÊNG LẺ (SO VỚI 0)
    for(int i=0;i<b;i++){
        cout <<"Điều kiện cho x"<<i+1<<" (1: ≥, 2: ≤, 3: tự do): ";
        while(true){
            cin >> dieukien_bien[i];
            if (dieukien_bien[i]>=1 && dieukien_bien[i]<=3) break;
            cout<<"Lỗi! Vui lòng nhập điều kiện hợp lệ (1, 2 hoặc 3): ";
        }
    }
    cout<<endl;
    
    int phphap;
    // CHỌN PHƯƠNG PHÁP GIẢI BÀI TOÁN QHTT
    cout<<"Chọn Phương pháp giải bài toán QHTT (1: Đơn hình, 2: Bland, 3: 2 pha, 4: Tất cả): ";
    while(true){
        cin>>phphap;
        if(phphap>=1 && phphap<=4)break;
        cout<<"Lỗi! Vui lòng nhập đúng phương pháp bạn muốn giải (1,2,3 hoặc 4): ";
    }
    
    
    
    // HÀM HIỂN THỊ MÔ HÌNH TỔNG QUÁT
    mohinh_tongquat(a,b,c,a_heso,rangbuoc_heso,rangbuoc_dau,sothucbi,dieukien_bien);
    cout<<"\n         ↓\n";
    // CHUYỂN ĐỔI SANG MÔ HÌNH CHUẨN TẮC VÀ LƯU TRỮ KẾT QUẢ
    // mo_hinh là đối tượng với kiểu MoHinhChuanTac
    MoHinhChuanTac mo_hinh = chuyensang_mohinhchuantac(a,b,c,a_heso, rangbuoc_heso, rangbuoc_dau, sothucbi, dieukien_bien);
    // HÀM HIỂN THỊ MÔ HÌNH CHUẨN TẮC SAU KHI ĐÃ CHUYỂN QUA TỪ MÔ HÌNH TỔNG QUÁT
    hienthi_MoHinhChuanTac(mo_hinh);
    // HIỂN THỊ TỪNG THÔNG TIN CỦA MÔ HÌNH CHUẨN TẮC
    mo_hinh.in_du_lieu();
    
    
    xu_ly_luachon(mo_hinh,phphap);
    
    
    // GIẢI PHÓNG DỮ LIỆU !!!
    delete[] a_heso;
    delete[] rangbuoc_dau;
    delete[] sothucbi;
    delete[] dieukien_bien;
    for (int i=0;i<c;i++){
        delete[] rangbuoc_heso[i];
    }
    delete[] rangbuoc_heso;

    return 0;
}