import os
import sys
import threading
import webbrowser
from src.app import app
import src.callbacks  # 必須載入以註冊 callbacks

def main():
    # 決定資源路徑 (支援 PyInstaller 打包後的路徑)
    if getattr(sys, 'frozen', False):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(os.path.abspath(__file__))

    # 自動開啟瀏覽器
    port = 8050
    threading.Timer(1.5, lambda: webbrowser.open(f"http://127.0.0.1:{port}")).start()
    
    print(f"✅ Pareto Explorer 啟動中...")
    print(f"🔗 請訪問 http://127.0.0.1:{port}")
    
    # 執行 App
    app.run(debug=False, port=port)

if __name__ == "__main__":
    main()
