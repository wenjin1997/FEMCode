- [引用参考文献](#引用参考文献)
  - [LaTeX在参考文献中引用网页](#latex在参考文献中引用网页)
- [列表使用](#列表使用)
- [插入Matlab代码](#插入matlab代码)
## 引用参考文献
* 实例：[数值分析笔记](https://github.com/wenjin1997/Numerical_Analysis)
* 具体见[参考文献使用文档](/biblatex-gb7714-2015.pdf)，[官方文档](https://github.com/hushidong/biblatex-gb7714-2015#jumptotutorial)
* 在自己的`document.tex`文档中加入

```
% 参考文献
% 注意参考文献请用Biber编译，不能使用BiberTex或者BiberTex-8
\usepackage[backend=biber,style=gb7714-2015ay]{biblatex}
\addbibresource[location=local]{bib/sample.bib}
```
* 文档开头引用包[可选]
```
%简单方式:
\usepackage[backend=biber,style=gb7714-2015]{biblatex} 
%设置gbalign选项以改变文献表序号标签对齐方式，设置gbpub=false取消缺省出版项自填补信息，比如: 
\usepackage[backend=biber,style=gb7714-2015,gbalign=gb7714-2015,gbpub=false]{biblatex} 
%当文档为GBK编码且用pdflatex/latex编译时，应设置选项gbcodegbk=true: 
\usepackage[backend=biber,style=gb7714-2015,gbcodegbk=true]{biblatex}
```
* 在`document.tex`文档末尾加上
```
\renewcommand\refname{参考文献}
%\bibliomatter
\nocite{*} %打印全部参考文献
\printbibliography[heading=bibliography,title=参考文献]
```
* 新建文件`/bib/sample.bib`，文档内容是参考文献，例如
```
%# -*- coding: utf-8-unix -*-

@Book{numerical_1,
  Title                    = {现代数值计算方法},
  Address                  = {北京},
  Author                   = {马昌凤 and 林伟川},
  Publisher                = {科学出版社},
  Year                     = {2008},
  Month                    = {6}
}

@Book{numerical_2,
  Title                    = {数值分析},
  Address                  = {北京},
  Author                   = {钟尔杰 and 黄廷祝},
  Publisher                = {高等教育出版社},
  Year                     = {2004},
  Month                    = {7}
}
```
* 编译顺序
```
xelatex document.tex
biber document
xelatex document.tex
xelatex document.tex
```
* 在正文中引用
```
\cite{numerical_1}
\citet{numerical_2}
```

### LaTeX在参考文献中引用网页
```
@misc{RN16,
  author = {Zanettin, Federico},
  year = {2000},
  url = {https://www.researchgate.net/publication/243771074_DIY_Corpora_the_WWW_and_the_Translator},
  urldate = {March 9, 2017},
  title = {DIY Corpora: the WWW and the Translator}
}
```
再引入url包
```
\usepackage{url}
```

## 列表使用
参考[LaTeX中列表的使用](https://blog.xulihang.me/use-list-in-latex/)
* 无序列表
```
\begin{itemize}
  \item This is the first item
  \item This is the second item
  \item This is the third item
\end{itemize}
```
* 有序列表
```
\begin{enumerate}
 \item This is the first item
 \item This is the second item
 \item This is the third item
\end{enumerate}
```

## 插入Matlab代码
参考[LaTeX中插入Matlab代码](https://blog.csdn.net/u012675539/article/details/47048163)