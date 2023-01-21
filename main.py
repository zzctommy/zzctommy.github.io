from flask import Flask, render_template, redirect, request
app = Flask(__name__)

data_chat = []
ip_user = {}
pre_message = {}

@app.route('/', methods = ["GET", "POST"])

def home():
    ip = request.remote_addr
    if not ip in ip_user:
        return redirect("/login")
    a = request.form.get("text")
    if a:
        if (ip in pre_message and a != pre_message[ip]) or not ip in pre_message:
            data_chat.append([ip, ip_user[ip], a])
        pre_message[ip] = a
        
    res = ""
    if len(data_chat) > 20:
        data_chat.pop(0)
    for i in data_chat:
        res += i[0] + ',' + i[1] + ' : '+ i[2] + '\n'
    return render_template('home.html', name = ip_user[ip], string = res)

@app.route("/login", methods = ["GET", "POST"])

def login():
    global is_login, username
    name = request.form.get("username")
    pwd = request.form.get("password")
    ip = request.remote_addr
    if name and pwd:
        username = name
        is_login = True
        ip_user[ip] = username
        return redirect('/')
    return render_template('login.html')

if __name__ == "__main__":
    app.run(host = "192.168.1.6", port = 5000)
