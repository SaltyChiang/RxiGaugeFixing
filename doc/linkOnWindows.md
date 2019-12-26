# PowerShell 创建软链接

    New-Item -ItemType SymbolicLink -Path path\\to\\link -Target path\\to\\source

    New-Item -ItemType HardLink -Path path\\to\\link -Target path\\to\\source

    New-Item -ItemType Junction -Path path\\to\\link -Target path\\to\\source

从上到下分别是软链接、硬链接和文件夹链接。