# bachelor_thesis_reproduce
中尾の卒論中の数値実験を再現するためのプログラムとデータ。


## sec3
論文第3章 "曲線座標系の要素を用いた有限要素法による解の数値検証"を再現する。


`model`ディレクトリおよび`output`ディレクトリでは、
`div1`,`div2`,`div3`,`div4`,`div5`がそれぞれ要素数
24,96,384,1536,6144のメッシュに対応する。

プログラムのコンパイル、実行は`sec3`ディレクトリで
```
make main && ./main
```
`visualize.ipynb`では、論文中図4の再現を確認する。

## sec4
論文第4章 "1次元問題への適用例"を再現する。


プログラムのコンパイル、実行は`sec4`ディレクトリで
```
make main && ./main
```
`visualize.ipynb`では、論文中図6,7,9,10,13,14の再現を確認する。

## sec5
論文第5章 "2次元問題への適用例"を再現する。

`ref/displacement_dist_cart_ref.csv`と`ref/stress_dist_cart_ref.csv`
はそれぞれ変位と応力の参照解であり、要素分割を細かくした2次要素を用いて求めたものである。

座標変換のパラメータ推定は、`sec5`ディレクトリで
```
python3 src/optimize_alpha.py
```
デカルト座標系での求めた解`output/displacement_node_cart.csv`が必要であることに注意。


また、プログラムのコンパイル、実行は`sec5`ディレクトリで
```
make main && ./main
```

`visualize.ipynb`では、論文中図16,17,18,19,20,21,22,23,24,25,26の再現を確認する。

## requirements
- gccコンパイラ(c++14対応)
- make
- python3インタプリタ
- python ライブラリ
  - jupyter
  - scipy
  - numpy
  - matplotlib
  - japanize_matplotlib
