# didfamily — 差の差(DiD)（jamovi モジュール）

同じパネルデータに対して複数の最新 DiD 推定量を並べて実行し、一致する点と食い違う点を確認します。診断、感度分析、再現可能な R コードも出力します。 詳しい英語の説明は [README.md](README.md)、中国語の要約は [README.zh.md](README.zh.md) を参照してください。

## 分析（jamovi メニュー：DiD）

| 分析 | 内容 |
|---|---|
| **パネル構造チェック** | 推定前にパネルを点検（均衡性、重複、欠損、処置の吸収性、未処置群の割合、コホートの大きさ） |
| **古典的 DiD (2x2)** | 2 群 2 期。群平均、回帰 DiD、二重頑健・IPW・結果回帰推定量（`DRDID`） |
| **推定量（個別分析）** | TWFE、Callaway & Sant'Anna、Sun & Abraham、Borusyak-Jaravel-Spiess（補完法）、Gardner 2 段階、拡張 TWFE（Wooldridge。線形・Poisson・logit）、積み上げ DiD。それぞれ全体 ATT、動的効果、イベントスタディ図 |
| **DiD 推定量の比較** | 最大 6 つの推定量を同じ交差処置パネルで比較。重ね描きのイベントスタディ図、フォレスト図、全数値を再現する R スクリプト |
| **Bacon 分解・事前トレンド・HonestDiD** | Goodman-Bacon 分解、TWFE ウェイト、事前トレンド検定と検出力、HonestDiD |
| **連続処置 DiD** | 用量反応 ATT(d) と限界効果 ACRT(d)（`contdid`） |
| **三重差分 (DDD)** | `triplediff` |
| **処置の ON/OFF 切り替え** | de Chaisemartin と D'Haultfoeuille (2024) の動的効果とプラセボ（独自実装、`DIDmultiplegtDYN` と一致） |

各分析の出力には「方法の説明」と「データと設定ガイド」があります。

標準誤差は既定で個体レベルのクラスター。約 300 のチェック（既知の効果を持つシミュレーション、`fixest`・`etwfe`・`bacondecomp` などとの一致）で検証しています。サンプルデータ（`.omv`）：`sim_staggered`、`sim_2x2`、`sim_ddd`、`sim_dose`、`sim_dcdh`、`sim_etwfe`、`mpdta`、`castle`。

## インストール
- **jamovi ライブラリ**：jamovi 右上のモジュールボタンから *jamovi library* を開き、検索してインストール（モジュール公開後）。
- **サイドロード**：`.jmo` ファイルをダウンロードし、同じメニューの *Install from file* から導入します。`.jmo` はビルドした OS、プロセッサ、jamovi のバージョンでのみ動作します。
- **ソースから**：R と `jmvtools` が必要です：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（パスはお使いの jamovi に合わせてください）。

## 引用とライセンス
引用：Sun, M. *didfamily* (jamovi module). https://github.com/mingchao-sun/didfamily 。使用した手法の原論文（出力に記載）も併せて引用してください。

ライセンス：`DESCRIPTION` を参照（GPL）。第三者のコードとデータは `NOTICE` に記載。不具合や提案：https://github.com/mingchao-sun/didfamily/issues
