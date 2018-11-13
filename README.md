# psico
## インストール方法 ##
まず `git clone https://github.com/BILAB/psico.git` などで中のpsicoディレクトリをダウンロードしたあと、これについてのPATHをPyMOL側から通しておきます。例えば`~/apps`にpsicoディレクトリを置いた場合は、ホームディレクトリ以下に`~/.pymolrc`を作成して
```
# ~/appsはpsicoディレクトリの場所に合わせる
sys.path.append(os.path.expanduser('~/apps'))
import psico.fullinit
```
としておきます。この状態でPyMOLを立ち上げれば、`contactmap`を始めとしたpsicoの中のプラグインに由来する追加コマンドが使えるようになっているはずです。
