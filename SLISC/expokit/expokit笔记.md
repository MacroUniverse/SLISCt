## ZGEXPV
Z 代表 Comp, G 代表 general matrix (不需要是 Hermitian)

```c++
void ZGEXPV(Int_I n, Int_I m, Doub_I t, CompVec_I &v, CompVec_O &w, Doub_I tol, Doub_I anorm, VecComp_O &wsp, Int_I lwsp, VecInt_O &iwsp, Int_I liwsp, FUNC matvec, Int_I itrace, Int_O &iflag)
```
其中 sparse matrix 的数据完全在函数 `matvec` 中, `ZGEXPV` 不需要知道矩阵的数据结构, 只需要通过 `matvec` 计算矩阵-矢量乘法.

* `n` 是方阵的尺寸
* `m` 是 Krylov basis 的数量 (?)
* `t` 是 `exp(At)` 中的 `t` (可以是负值)
* `v(n)` 和 `w(n)` 是 `w=exp(At)v` 中的.
* `tol` 是 `w` 的精度， 最小值是 `sqrt(eps)` (输入 `tol=0`).
* `anorm` 用于衡量矩阵元数值的大小, 例程中用的是 `anorm = max(sum(abs(a),2))`
* `wsp(lwsp)` 是 Comp 类型的 work space, `lwsp >= n*(m+2)+5*(m+2)^2+ideg+1` 其中 `ideg=6` 是程序的内部常数.
* `iwsp(liwsp)` 是 Int 类型的 work space, `liwsp >= m+2`
* `itrace` 输出控制, 0=silent, 1=print step-by-step info
* `iflag` exit flag. `<0`: bad input arguments, `0`: no problem, `1`: maximum number of steps reached without convergence, `2`: requested tolerance was too high.
* 程序执行完以后, `wsp` 和 `iwsp` 里面会有程序运行的信息, 见程序的注释. 如果不需要可以从程序里面删除以提升速度.

## ZHEXPV
与 `ZGEXPV` 格式一样, 只是 `H` 代表 Hermitian 矩阵.


## ZGCOOV ZGCRSV ZGCCSV
计算 COO, CRS, CCS 形式的的 sparse matrix 和列矢量的乘法, 
```c++
void ZG???V(VecComp_I &x, VecComp_O &y)
```
`ZGCOOV` 中 sparse matrix 通过 `/common/ a, ia, ja, nz, n` 来传递, `nz` 非零矩阵元的个数, `n` 是方阵的尺寸.
* 如果写成 c++, 可以给 COO/CRS/CCS 做一个 class, 然后定义乘法 method 即可. 这样就可以通过 reference 传递给 `ZGEXPV` 等函数.
