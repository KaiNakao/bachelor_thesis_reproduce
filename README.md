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