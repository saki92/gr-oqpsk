git add .
git reset -- build/
echo "Enter a commit message:"
read cmt_msg
git commit -m "$cmt_msg"
